//! Shared topology projection used by object-model adapters.

use pdbiox_chem::element_properties;
use pdbiox_core::contract::Namespace;
use pdbiox_core::structure::{AtomRef, ChainRef, ResidueRef};
use pdbiox_core::{BondOrder, ModelIndex, Structure};
use thiserror::Error;

/// A chain in an external object-model projection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExportChain {
    /// Normalised chain identifier.
    pub id: String,
    /// Half-open residue range.
    pub residues: std::ops::Range<usize>,
}

/// A residue in an external object-model projection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExportResidue {
    /// Component identifier.
    pub name: String,
    /// Depositor residue number, when present.
    pub number: Option<i32>,
    /// Insertion code, when present.
    pub insertion_code: Option<String>,
    /// Whether the source marked this residue as a heterogen.
    pub is_heterogen: bool,
    /// Parent chain position.
    pub chain: usize,
    /// Half-open atom range.
    pub atoms: std::ops::Range<usize>,
}

/// An atom in an external object-model projection.
#[derive(Clone, Debug, PartialEq)]
pub struct ExportAtom {
    /// Atom name.
    pub name: String,
    /// Atomic number.
    pub atomic_number: u8,
    /// Canonical element symbol.
    pub element_symbol: String,
    /// Standard atomic weight from the chemistry table.
    pub mass: f64,
    /// File-local atom identifier, when present.
    pub serial: Option<u32>,
    /// Formal charge, when present.
    pub formal_charge: Option<i8>,
    /// Occupancy, when present.
    pub occupancy: Option<f32>,
    /// Temperature factor, when present.
    pub b_factor: Option<f32>,
    /// Alternate-location identifier, when present.
    pub alternate_location: Option<String>,
    /// Parent residue position.
    pub residue: usize,
    /// Cartesian position in ångström.
    pub position: [f32; 3],
}

/// A bond in an external object-model projection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExportBond {
    /// First atom position.
    pub atom_a: usize,
    /// Second atom position.
    pub atom_b: usize,
    /// Chemical order.
    pub order: BondOrder,
}

/// A complete, dense projection produced in one native pass.
#[derive(Clone, Debug, PartialEq)]
pub struct TopologyExport {
    /// Chains in structure order.
    pub chains: Vec<ExportChain>,
    /// Residues in structure order.
    pub residues: Vec<ExportResidue>,
    /// Atoms in structure order.
    pub atoms: Vec<ExportAtom>,
    /// Bonds in endpoint order.
    pub bonds: Vec<ExportBond>,
}

/// Why an object-model projection cannot be made without inventing data.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum TopologyExportError {
    /// Explicit namespaces require a caller-owned mapping unavailable here.
    #[error("topology export requires label or auth identifiers")]
    ExplicitNamespace,
    /// The requested model does not exist.
    #[error("model {model} does not exist")]
    ModelUnavailable {
        /// Requested zero-based model position.
        model: usize,
    },
    /// A required structural field is absent.
    #[error("{kind} {index} has no {field}")]
    MissingField {
        /// Entity kind.
        kind: &'static str,
        /// Entity position.
        index: usize,
        /// Missing field.
        field: &'static str,
    },
    /// A hierarchy range does not tile the projected rows.
    #[error("{kind} {index} is not contiguous in the structure hierarchy")]
    NonContiguousHierarchy {
        /// Entity kind.
        kind: &'static str,
        /// Entity position.
        index: usize,
    },
    /// Connectivity refers outside the projected atoms.
    #[error("bond endpoint {atom} is outside {atom_count} projected atoms")]
    InvalidBondEndpoint {
        /// Invalid atom position.
        atom: usize,
        /// Number of projected atoms.
        atom_count: usize,
    },
}

impl TopologyExport {
    /// Projects one model without guessing missing structural data.
    ///
    /// # Errors
    ///
    /// Returns an explicit error when the model is unavailable, required
    /// topology or chemistry is absent, or hierarchy/connectivity is invalid.
    pub fn from_model(
        structure: &Structure,
        model: ModelIndex,
        namespace: Namespace,
    ) -> Result<Self, TopologyExportError> {
        if namespace == Namespace::Explicit {
            return Err(TopologyExportError::ExplicitNamespace);
        }
        let requested = model.as_usize();
        let (snapshot, local_model) = structure
            .model_snapshot(model)
            .ok_or(TopologyExportError::ModelUnavailable { model: requested })?;
        let positions = snapshot
            .model_positions(local_model)
            .ok_or(TopologyExportError::ModelUnavailable { model: requested })?;
        let data = snapshot.data();
        let mut export = Self {
            chains: Vec::with_capacity(snapshot.chain_count()),
            residues: Vec::with_capacity(snapshot.residue_count()),
            atoms: Vec::with_capacity(snapshot.atom_count() as usize),
            bonds: Vec::with_capacity(data.bonds.len()),
        };
        export.project_hierarchy(&snapshot, positions, namespace)?;
        export.project_bonds(&snapshot)?;
        Ok(export)
    }

    fn project_hierarchy(
        &mut self,
        structure: &Structure,
        positions: &[[f32; 3]],
        namespace: Namespace,
    ) -> Result<(), TopologyExportError> {
        for (chain_position, chain) in structure.data().chains().enumerate() {
            let residue_start = self.residues.len();
            let id = required(
                chain_id(chain, namespace),
                "chain",
                chain_position,
                "identifier",
            )?;
            for residue in chain.residues() {
                let residue_position = self.residues.len();
                let atom_start = self.atoms.len();
                let name = required(
                    residue_name(residue, namespace),
                    "residue",
                    residue_position,
                    "name",
                )?;
                for atom in residue.atoms() {
                    let atom_position = self.atoms.len();
                    if atom.index().as_usize() != atom_position {
                        return Err(TopologyExportError::NonContiguousHierarchy {
                            kind: "atom",
                            index: atom_position,
                        });
                    }
                    let element = required(atom.element(), "atom", atom_position, "element")?;
                    let properties = required(
                        element_properties(element),
                        "atom",
                        atom_position,
                        "element properties",
                    )?;
                    let position = required(
                        positions.get(atom_position).copied(),
                        "atom",
                        atom_position,
                        "coordinates",
                    )?;
                    self.atoms.push(ExportAtom {
                        name: required(atom_name(atom, namespace), "atom", atom_position, "name")?
                            .to_owned(),
                        atomic_number: element.atomic_number(),
                        element_symbol: element.symbol().to_owned(),
                        mass: properties.atomic_weight,
                        serial: atom.atom_site_id(),
                        formal_charge: atom.formal_charge(),
                        occupancy: atom.occupancy(),
                        b_factor: atom.b_factor(),
                        alternate_location: atom.alt_label().map(str::to_owned),
                        residue: residue_position,
                        position,
                    });
                }
                self.residues.push(ExportResidue {
                    name: name.to_owned(),
                    number: residue_number(residue, namespace),
                    insertion_code: residue.ins_code().map(str::to_owned),
                    is_heterogen: residue.is_het(),
                    chain: chain_position,
                    atoms: atom_start..self.atoms.len(),
                });
            }
            self.chains.push(ExportChain {
                id: id.to_owned(),
                residues: residue_start..self.residues.len(),
            });
        }
        if self.atoms.len() != positions.len() {
            return Err(TopologyExportError::NonContiguousHierarchy {
                kind: "atom",
                index: self.atoms.len(),
            });
        }
        Ok(())
    }

    fn project_bonds(&mut self, structure: &Structure) -> Result<(), TopologyExportError> {
        for bond in structure.data().bonds.iter() {
            let atom_a = bond.atom_a.as_usize();
            let atom_b = bond.atom_b.as_usize();
            for endpoint in [atom_a, atom_b] {
                if endpoint >= self.atoms.len() {
                    return Err(TopologyExportError::InvalidBondEndpoint {
                        atom: endpoint,
                        atom_count: self.atoms.len(),
                    });
                }
            }
            self.bonds.push(ExportBond {
                atom_a,
                atom_b,
                order: bond.order,
            });
        }
        Ok(())
    }
}

fn chain_id(chain: ChainRef<'_>, namespace: Namespace) -> Option<&str> {
    match namespace {
        Namespace::Label => chain.label(),
        Namespace::Auth => chain.auth_label(),
        _ => None,
    }
}

fn residue_name(residue: ResidueRef<'_>, namespace: Namespace) -> Option<&str> {
    match namespace {
        Namespace::Label => residue.name(),
        Namespace::Auth => residue.auth_name(),
        _ => None,
    }
}

fn residue_number(residue: ResidueRef<'_>, namespace: Namespace) -> Option<i32> {
    match namespace {
        Namespace::Label => residue.label_seq_id(),
        Namespace::Auth => residue.auth_seq_id(),
        _ => None,
    }
}

fn atom_name(atom: AtomRef<'_>, namespace: Namespace) -> Option<&str> {
    match namespace {
        Namespace::Label => atom.name(),
        Namespace::Auth => atom.auth_name(),
        _ => None,
    }
}

fn required<T>(
    value: Option<T>,
    kind: &'static str,
    index: usize,
    field: &'static str,
) -> Result<T, TopologyExportError> {
    value.ok_or(TopologyExportError::MissingField { kind, index, field })
}

#[cfg(test)]
#[path = "topology_tests.rs"]
mod tests;
