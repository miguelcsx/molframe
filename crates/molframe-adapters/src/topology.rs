//! Compact columnar topology projection for neutral interoperability.

use molframe_chem::element_properties;
use molframe_core::column::{BitVec, Presence, ValidityMask};
use molframe_core::contract::Namespace;
use molframe_core::structure::{AtomRef, ChainRef, ResidueRef};
use molframe_core::{BondOrder, ModelIndex, Structure};
use std::collections::HashMap;
use thiserror::Error;

/// Sentinel for an absent entry in the batch-local string dictionary.
pub const MISSING_STRING: u32 = u32::MAX;

/// A topology snapshot stored entirely as contiguous columns.
#[derive(Clone, Debug, PartialEq)]
pub struct TopologyBatch {
    /// Shared strings referenced by compact column IDs.
    pub strings: Vec<String>,
    /// Chain identifier dictionary IDs.
    pub chain_ids: Vec<u32>,
    /// Chain-to-residue offsets; length is `chain_count + 1`.
    pub chain_residue_offsets: Vec<u32>,
    /// Residue name dictionary IDs.
    pub residue_names: Vec<u32>,
    /// Depositor residue numbers; consult `residue_number_validity`.
    pub residue_numbers: Vec<i32>,
    /// Validity for `residue_numbers`.
    pub residue_number_validity: ValidityMask,
    /// Insertion-code dictionary IDs or [`MISSING_STRING`].
    pub residue_insertion_codes: Vec<u32>,
    /// Compact heterogen flags.
    pub residue_is_heterogen: BitVec,
    /// Parent chain per residue.
    pub residue_chain: Vec<u32>,
    /// Residue-to-atom offsets; length is `residue_count + 1`.
    pub residue_atom_offsets: Vec<u32>,
    /// Atom-name dictionary IDs.
    pub atom_names: Vec<u32>,
    /// Atomic numbers.
    pub atomic_numbers: Vec<u8>,
    /// Standard atomic weights.
    pub masses: Vec<f64>,
    /// File-local serials; consult `atom_serial_validity`.
    pub atom_serials: Vec<u32>,
    /// Validity for `atom_serials`.
    pub atom_serial_validity: ValidityMask,
    /// Formal charges; consult `formal_charge_validity`.
    pub formal_charges: Vec<i8>,
    /// Validity for `formal_charges`.
    pub formal_charge_validity: ValidityMask,
    /// Occupancies; consult `occupancy_validity`.
    pub occupancies: Vec<f32>,
    /// Validity for `occupancies`.
    pub occupancy_validity: ValidityMask,
    /// Temperature factors; consult `b_factor_validity`.
    pub b_factors: Vec<f32>,
    /// Validity for `b_factors`.
    pub b_factor_validity: ValidityMask,
    /// Alternate-location dictionary IDs or [`MISSING_STRING`].
    pub atom_alternate_locations: Vec<u32>,
    /// Parent residue per atom.
    pub atom_residue: Vec<u32>,
    /// Cartesian x coordinates.
    pub position_x: Vec<f32>,
    /// Cartesian y coordinates.
    pub position_y: Vec<f32>,
    /// Cartesian z coordinates.
    pub position_z: Vec<f32>,
    /// First bond endpoints.
    pub bond_atom_a: Vec<u32>,
    /// Second bond endpoints.
    pub bond_atom_b: Vec<u32>,
    /// Bond orders aligned with endpoints.
    pub bond_orders: Vec<BondOrder>,
}

/// Why a topology batch cannot be produced without inventing data.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum TopologyBatchError {
    /// Explicit namespaces require a caller-owned mapping unavailable here.
    #[error("topology export requires label or auth identifiers")]
    ExplicitNamespace,
    /// The requested model does not exist.
    #[error("model {model} does not exist")]
    ModelUnavailable {
        /// Requested zero-based model.
        model: usize,
    },
    /// A required structural field is absent.
    #[error("{kind} {index} has no {field}")]
    MissingField {
        /// Row family.
        kind: &'static str,
        /// Row position.
        index: usize,
        /// Required field.
        field: &'static str,
    },
    /// A hierarchy range does not tile the projected rows.
    #[error("{kind} {index} is not contiguous in the structure hierarchy")]
    NonContiguousHierarchy {
        /// Row family.
        kind: &'static str,
        /// First non-contiguous row.
        index: usize,
    },
    /// A row count cannot fit the compact local index space.
    #[error("{kind} count {count} exceeds the local row space")]
    RowCapacity {
        /// Row family.
        kind: &'static str,
        /// Requested rows.
        count: usize,
    },
    /// Connectivity refers outside the projected atoms.
    #[error("bond endpoint {atom} is outside {atom_count} projected atoms")]
    InvalidBondEndpoint {
        /// Invalid atom position.
        atom: usize,
        /// Available atoms.
        atom_count: usize,
    },
}

impl TopologyBatch {
    /// Projects one bounded model into compact columns in one pass.
    ///
    /// # Errors
    ///
    /// Returns an explicit error for an unavailable model, missing required
    /// topology, invalid connectivity, or local index overflow.
    pub fn from_model(
        structure: &Structure,
        model: ModelIndex,
        namespace: Namespace,
    ) -> Result<Self, TopologyBatchError> {
        if namespace == Namespace::Explicit {
            return Err(TopologyBatchError::ExplicitNamespace);
        }
        let requested = model.as_usize();
        let (snapshot, local_model) = structure
            .model_snapshot(model)
            .ok_or(TopologyBatchError::ModelUnavailable { model: requested })?;
        let positions = snapshot
            .model_positions(local_model)
            .ok_or(TopologyBatchError::ModelUnavailable { model: requested })?;
        let atom_count = snapshot.atom_count() as usize;
        let residue_count = snapshot.residue_count();
        let mut output = Self::with_capacity(
            snapshot.chain_count(),
            residue_count,
            atom_count,
            snapshot.data().bonds.len(),
        )?;
        let mut dictionary = HashMap::new();
        output.project_hierarchy(&snapshot, positions, namespace, &mut dictionary)?;
        output.project_bonds(&snapshot)?;
        Ok(output)
    }

    /// Resolves a batch-local string ID.
    #[must_use]
    pub fn string(&self, id: u32) -> Option<&str> {
        self.strings.get(id as usize).map(String::as_str)
    }

    fn with_capacity(
        chains: usize,
        residues: usize,
        atoms: usize,
        bonds: usize,
    ) -> Result<Self, TopologyBatchError> {
        let atom_len = local_len("atom", atoms)?;
        let residue_len = local_len("residue", residues)?;
        let _chain_len = local_len("chain", chains)?;
        Ok(Self {
            strings: Vec::new(),
            chain_ids: Vec::with_capacity(chains),
            chain_residue_offsets: vec![0],
            residue_names: Vec::with_capacity(residues),
            residue_numbers: Vec::with_capacity(residues),
            residue_number_validity: ValidityMask::all_present(residue_len),
            residue_insertion_codes: Vec::with_capacity(residues),
            residue_is_heterogen: BitVec::repeat(false, residue_len),
            residue_chain: Vec::with_capacity(residues),
            residue_atom_offsets: vec![0],
            atom_names: Vec::with_capacity(atoms),
            atomic_numbers: Vec::with_capacity(atoms),
            masses: Vec::with_capacity(atoms),
            atom_serials: Vec::with_capacity(atoms),
            atom_serial_validity: ValidityMask::all_present(atom_len),
            formal_charges: Vec::with_capacity(atoms),
            formal_charge_validity: ValidityMask::all_present(atom_len),
            occupancies: Vec::with_capacity(atoms),
            occupancy_validity: ValidityMask::all_present(atom_len),
            b_factors: Vec::with_capacity(atoms),
            b_factor_validity: ValidityMask::all_present(atom_len),
            atom_alternate_locations: Vec::with_capacity(atoms),
            atom_residue: Vec::with_capacity(atoms),
            position_x: Vec::with_capacity(atoms),
            position_y: Vec::with_capacity(atoms),
            position_z: Vec::with_capacity(atoms),
            bond_atom_a: Vec::with_capacity(bonds),
            bond_atom_b: Vec::with_capacity(bonds),
            bond_orders: Vec::with_capacity(bonds),
        })
    }

    fn project_hierarchy(
        &mut self,
        structure: &Structure,
        positions: &[[f32; 3]],
        namespace: Namespace,
        dictionary: &mut HashMap<String, u32>,
    ) -> Result<(), TopologyBatchError> {
        for (chain_position, chain) in structure.data().chains().enumerate() {
            let id = required(
                chain_id(chain, namespace),
                "chain",
                chain_position,
                "identifier",
            )?;
            self.chain_ids
                .push(intern(&mut self.strings, dictionary, id)?);
            for residue in chain.residues() {
                self.push_residue(residue, chain_position, namespace, dictionary)?;
                for atom in residue.atoms() {
                    self.push_atom(atom, positions, namespace, dictionary)?;
                }
                self.residue_atom_offsets
                    .push(local_len("atom", self.atom_names.len())?);
            }
            self.chain_residue_offsets
                .push(local_len("residue", self.residue_names.len())?);
        }
        if self.atom_names.len() != positions.len() {
            return Err(TopologyBatchError::NonContiguousHierarchy {
                kind: "atom",
                index: self.atom_names.len(),
            });
        }
        Ok(())
    }

    fn push_residue(
        &mut self,
        residue: ResidueRef<'_>,
        chain: usize,
        namespace: Namespace,
        dictionary: &mut HashMap<String, u32>,
    ) -> Result<(), TopologyBatchError> {
        let index = self.residue_names.len();
        let name = required(residue_name(residue, namespace), "residue", index, "name")?;
        self.residue_names
            .push(intern(&mut self.strings, dictionary, name)?);
        push_optional(
            &mut self.residue_numbers,
            &mut self.residue_number_validity,
            index,
            residue_number(residue, namespace),
        );
        self.residue_insertion_codes.push(match residue.ins_code() {
            Some(value) => intern(&mut self.strings, dictionary, value)?,
            None => MISSING_STRING,
        });
        self.residue_is_heterogen
            .set(local_len("residue", index)?, residue.is_het());
        self.residue_chain.push(local_len("chain", chain)?);
        Ok(())
    }

    fn push_atom(
        &mut self,
        atom: AtomRef<'_>,
        positions: &[[f32; 3]],
        namespace: Namespace,
        dictionary: &mut HashMap<String, u32>,
    ) -> Result<(), TopologyBatchError> {
        let index = self.atom_names.len();
        if atom.index().as_usize() != index {
            return Err(TopologyBatchError::NonContiguousHierarchy {
                kind: "atom",
                index,
            });
        }
        let element = required(atom.element(), "atom", index, "element")?;
        let properties = required(
            element_properties(element),
            "atom",
            index,
            "element properties",
        )?;
        let position = required(positions.get(index).copied(), "atom", index, "coordinates")?;
        let name = required(atom_name(atom, namespace), "atom", index, "name")?;
        self.atom_names
            .push(intern(&mut self.strings, dictionary, name)?);
        self.atomic_numbers.push(element.atomic_number());
        self.masses.push(properties.atomic_weight);
        push_optional(
            &mut self.atom_serials,
            &mut self.atom_serial_validity,
            index,
            atom.atom_site_id(),
        );
        push_optional(
            &mut self.formal_charges,
            &mut self.formal_charge_validity,
            index,
            atom.formal_charge(),
        );
        push_optional(
            &mut self.occupancies,
            &mut self.occupancy_validity,
            index,
            atom.occupancy(),
        );
        push_optional(
            &mut self.b_factors,
            &mut self.b_factor_validity,
            index,
            atom.b_factor(),
        );
        self.atom_alternate_locations.push(match atom.alt_label() {
            Some(value) => intern(&mut self.strings, dictionary, value)?,
            None => MISSING_STRING,
        });
        self.atom_residue
            .push(local_len("residue", self.residue_names.len() - 1)?);
        self.position_x.push(position[0]);
        self.position_y.push(position[1]);
        self.position_z.push(position[2]);
        Ok(())
    }

    fn project_bonds(&mut self, structure: &Structure) -> Result<(), TopologyBatchError> {
        for bond in structure.data().bonds.iter() {
            let atom_a = bond.atom_a.get();
            let atom_b = bond.atom_b.get();
            for endpoint in [atom_a, atom_b] {
                if endpoint as usize >= self.atom_names.len() {
                    return Err(TopologyBatchError::InvalidBondEndpoint {
                        atom: endpoint as usize,
                        atom_count: self.atom_names.len(),
                    });
                }
            }
            self.bond_atom_a.push(atom_a);
            self.bond_atom_b.push(atom_b);
            self.bond_orders.push(bond.order);
        }
        Ok(())
    }
}

fn intern(
    strings: &mut Vec<String>,
    dictionary: &mut HashMap<String, u32>,
    value: &str,
) -> Result<u32, TopologyBatchError> {
    if let Some(id) = dictionary.get(value) {
        return Ok(*id);
    }
    let id = local_len("string dictionary", strings.len())?;
    strings.push(value.to_owned());
    dictionary.insert(value.to_owned(), id);
    Ok(id)
}

fn push_optional<T: Copy + Default>(
    values: &mut Vec<T>,
    validity: &mut ValidityMask,
    index: usize,
    value: Option<T>,
) {
    if let Some(value) = value {
        values.push(value);
    } else {
        values.push(T::default());
        if let Ok(position) = u32::try_from(index) {
            validity.set(position, Presence::Inapplicable);
        }
    }
}

fn local_len(kind: &'static str, count: usize) -> Result<u32, TopologyBatchError> {
    u32::try_from(count).map_err(|_| TopologyBatchError::RowCapacity { kind, count })
}

fn chain_id(chain: ChainRef<'_>, namespace: Namespace) -> Option<&str> {
    match namespace {
        Namespace::Label => chain.label(),
        Namespace::Auth => chain.auth_label().or_else(|| chain.label()),
        _ => None,
    }
}

fn residue_name(residue: ResidueRef<'_>, namespace: Namespace) -> Option<&str> {
    match namespace {
        Namespace::Label => residue.name(),
        Namespace::Auth => residue.auth_name().or_else(|| residue.name()),
        _ => None,
    }
}

fn residue_number(residue: ResidueRef<'_>, namespace: Namespace) -> Option<i32> {
    match namespace {
        Namespace::Label => residue.label_seq_id(),
        Namespace::Auth => residue.auth_seq_id().or_else(|| residue.label_seq_id()),
        _ => None,
    }
}

fn atom_name(atom: AtomRef<'_>, namespace: Namespace) -> Option<&str> {
    match namespace {
        Namespace::Label => atom.name(),
        Namespace::Auth => atom.auth_name().or_else(|| atom.name()),
        _ => None,
    }
}

fn required<T>(
    value: Option<T>,
    kind: &'static str,
    index: usize,
    field: &'static str,
) -> Result<T, TopologyBatchError> {
    value.ok_or(TopologyBatchError::MissingField { kind, index, field })
}

#[cfg(test)]
#[path = "topology_tests.rs"]
mod tests;
