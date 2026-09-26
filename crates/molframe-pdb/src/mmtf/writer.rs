use super::codec::{encode_chars, encode_f32, encode_i32, encode_strings};
use super::schema::{Entity, File, Group};
use super::{MMTF_METADATA_EXTENSION, MmtfMetadata, MmtfOptionalField};
use molframe_core::BondOrder;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::{EntityIndex, ModelIndex};
use molframe_core::structure::Structure;
use num_traits::ToPrimitive;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::io::Write;

const VERSION: &str = "1.0.0";
const PRODUCER: &str = "molframe 0.1.1";
const CHAIN_WIDTH: usize = 4;

/// Writes a deterministic standards-conforming MMTF `MessagePack` document.
///
/// # Errors
///
/// Refuses ragged models, absent coordinates, identifiers beyond MMTF's byte
/// limits, non-ASCII character fields, and counts outside signed 32-bit range.
pub fn write_mmtf(structure: &Structure) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let mut output = Vec::new();
    write_mmtf_to(structure, &mut output)?;
    Ok(output)
}

/// Streams deterministic MMTF `MessagePack` to a byte destination.
///
/// # Errors
///
/// Returns structure, encoding, or destination diagnostics.
pub fn write_mmtf_to<W: Write>(
    structure: &Structure,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    validate_structure(structure).map_err(|finding| vec![finding])?;
    let file = build(structure).map_err(|finding| vec![finding])?;
    rmp_serde::encode::write_named(output, &file).map_err(|error| {
        vec![
            Diagnostic::new(Code::E1102)
                .with_message("MMTF MessagePack could not be encoded")
                .with_context("encoder", error.to_string()),
        ]
    })
}

fn validate_structure(structure: &Structure) -> Result<(), Diagnostic> {
    if structure.ragged_models().is_some() {
        Err(refusal(
            "MMTF writer requires shared topology across models",
        ))
    } else {
        Ok(())
    }
}

fn build(structure: &Structure) -> Result<File, Diagnostic> {
    let models = structure.model_count();
    let atom_count = structure.atom_count() as usize;
    let residue_count = structure.residue_count();
    let chain_count = structure.chain_count();
    let metadata = source_metadata(structure)?;
    let mut group_list = Vec::with_capacity(residue_count);
    let mut topology = TopologyColumns::with_capacity(residue_count, chain_count, models);
    let atom_capacity = atom_count
        .checked_mul(models)
        .ok_or_else(|| refusal("MMTF atom capacity overflows"))?;
    let mut atoms = AtomColumns::with_capacity(atom_capacity);
    let group_type_of = build_group_dictionary(structure, metadata, &mut group_list)?;
    append_models(structure, models, &group_type_of, &mut topology, &mut atoms)?;
    let (bond_atoms, bond_orders, bond_resonance) = bonds(structure, models)?;
    let num_bonds = i32_of(bond_orders.len(), "bond count")?;
    let (bond_atom_list, bond_order_list, bond_resonance_list) =
        encode_bond_columns(&bond_atoms, &bond_orders, &bond_resonance)?;
    Ok(File {
        mmtf_version: VERSION.to_string(),
        mmtf_producer: PRODUCER.to_string(),
        structure_id: structure.data().entry.id.as_deref().map(str::to_owned),
        title: structure.data().entry.title.as_deref().map(str::to_owned),
        unit_cell: encoded_cell(structure)?,
        space_group: metadata.space_group.as_deref().map(str::to_owned),
        experimental_methods: experimental_methods(structure),
        resolution: structure.data().entry.resolution,
        num_bonds,
        num_atoms: i32_of(atoms.x.len(), "atom count")?,
        num_groups: i32_of(topology.group_ids.len(), "group count")?,
        num_chains: i32_of(topology.chain_ids.len(), "chain count")?,
        num_models: i32_of(models, "model count")?,
        group_list,
        bond_atom_list,
        bond_order_list,
        bond_resonance_list,
        x_coord_list: encode_f32(&atoms.x)?,
        y_coord_list: encode_f32(&atoms.y)?,
        z_coord_list: encode_f32(&atoms.z)?,
        b_factor_list: encode_optional_f32(
            atoms.b_factor.values(),
            metadata
                .optional_fields
                .contains(&MmtfOptionalField::BFactor),
            "B-factor list",
        )?,
        atom_id_list: encode_optional_i32(
            atoms.atom_id.values(),
            metadata
                .optional_fields
                .contains(&MmtfOptionalField::AtomId),
            "atom identifier list",
        )?,
        alt_loc_list: chars_with_presence(
            &atoms.alt_loc,
            metadata
                .optional_fields
                .contains(&MmtfOptionalField::AltLoc),
            "alternate-location list",
        )?,
        occupancy_list: encode_optional_f32(
            atoms.occupancy.values(),
            metadata
                .optional_fields
                .contains(&MmtfOptionalField::Occupancy),
            "occupancy list",
        )?,
        group_id_list: encode_i32(&topology.group_ids)?,
        group_type_list: encode_i32(&topology.group_types)?,
        ins_code_list: chars_with_presence(
            &topology.insertion_codes,
            metadata
                .optional_fields
                .contains(&MmtfOptionalField::InsCode),
            "insertion-code list",
        )?,
        sequence_index_list: metadata
            .optional_fields
            .contains(&MmtfOptionalField::SequenceIndex)
            .then(|| encode_i32(&topology.sequence_indices))
            .transpose()?,
        chain_id_list: encode_strings(&topology.chain_ids, CHAIN_WIDTH)?,
        chain_name_list: required_presence(
            topology.chain_names.values(),
            metadata
                .optional_fields
                .contains(&MmtfOptionalField::ChainName),
            "author chain name list",
        )?
        .map(|names| encode_strings(names, CHAIN_WIDTH))
        .transpose()?,
        groups_per_chain: topology.groups_per_chain,
        chains_per_model: vec![i32_of(chain_count, "chains per model")?; models],
        entity_list: metadata
            .optional_fields
            .contains(&MmtfOptionalField::EntityList)
            .then(|| entities(structure, metadata, models))
            .transpose()?,
    })
}

fn append_models(
    structure: &Structure,
    models: usize,
    group_type_of: &[i32],
    topology: &mut TopologyColumns,
    atoms: &mut AtomColumns,
) -> Result<(), Diagnostic> {
    for model in 0..models {
        topology.append(structure, group_type_of)?;
        let model = u32::try_from(model).map_err(|_| refusal("MMTF model index exceeds u32"))?;
        atoms.append(structure, ModelIndex::new(model))?;
    }
    Ok(())
}

/// Builds the group dictionary and the residue-to-group-type mapping.
///
/// MMTF stores each distinct chemistry once and refers to it by index, so a
/// structure with a million glycines carries one glycine entry. Emitting one
/// entry per residue instead made the dictionary as large as the structure and
/// repeated every name and element list with it.
fn build_group_dictionary(
    structure: &Structure,
    metadata: &MmtfMetadata,
    groups: &mut Vec<Group>,
) -> Result<Vec<i32>, Diagnostic> {
    // The metadata chemistry list is searched once per residue, so it is keyed
    // by component name first rather than scanned end to end each time.
    let mut chemistry_by_name: HashMap<&str, Vec<usize>> = HashMap::new();
    for (position, candidate) in metadata.groups.iter().enumerate() {
        chemistry_by_name
            .entry(candidate.name.as_ref())
            .or_default()
            .push(position);
    }

    let mut group_type_of = vec![0i32; structure.residue_count()];
    let mut seen: HashMap<Group, i32> = HashMap::new();

    for residue in structure.data().residues() {
        let atoms: Vec<_> = residue.atoms().collect();
        let formal_charge_list = atoms
            .iter()
            .map(|atom| {
                atom.formal_charge()
                    .map(i32::from)
                    .ok_or_else(|| refusal("MMTF requires a formal charge for every atom"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let atom_name_list = atoms
            .iter()
            .map(|atom| required_text(atom.name(), "atom name"))
            .collect::<Result<Vec<_>, _>>()?;
        let group_name = required_text(residue.name(), "component name")?;
        let chemistry = chemistry_by_name
            .get(group_name.as_str())
            .and_then(|candidates| {
                candidates.iter().find_map(|position| {
                    let candidate = metadata.groups.get(*position)?;
                    (candidate.atom_names.len() == atom_name_list.len()
                        && candidate
                            .atom_names
                            .iter()
                            .zip(&atom_name_list)
                            .all(|(expected, actual)| expected.as_ref() == actual))
                    .then_some(candidate)
                })
            })
            .ok_or_else(|| refusal("MMTF group chemistry no longer matches the structure"))?;
        let element_list = chemistry
            .elements
            .as_ref()
            .map(|elements| validate_elements(&atoms, elements))
            .transpose()?;
        if atom_name_list.iter().any(|name| name.len() > 5) || group_name.len() > 5 {
            return Err(refusal(
                "MMTF atom and component names are limited to five bytes",
            ));
        }
        let group = Group {
            formal_charge_list,
            atom_name_list,
            element_list,
            bond_atom_list: Vec::new(),
            bond_order_list: Vec::new(),
            bond_resonance_list: Vec::new(),
            name: group_name,
            single_letter_code: chemistry.single_letter_code.to_string(),
            chem_comp_type: chemistry.chem_comp_type.to_string(),
        };

        // First appearance order decides the identifiers, so the dictionary is
        // the same for the same structure on every run.
        let next = i32_of(groups.len(), "group type")?;
        let group_type = match seen.entry(group) {
            Entry::Occupied(entry) => *entry.get(),
            Entry::Vacant(entry) => {
                groups.push(entry.key().clone());
                *entry.insert(next)
            }
        };

        if let Some(slot) = group_type_of.get_mut(residue.index().as_usize()) {
            *slot = group_type;
        }
    }
    Ok(group_type_of)
}

struct TopologyColumns {
    group_ids: Vec<i32>,
    group_types: Vec<i32>,
    insertion_codes: Vec<u8>,
    sequence_indices: Vec<i32>,
    groups_per_chain: Vec<i32>,
    chain_ids: Vec<String>,
    chain_names: OptionalColumn<String>,
}

impl TopologyColumns {
    fn with_capacity(residues: usize, chains: usize, models: usize) -> Self {
        Self {
            group_ids: Vec::with_capacity(residues * models),
            group_types: Vec::with_capacity(residues * models),
            insertion_codes: Vec::with_capacity(residues * models),
            sequence_indices: Vec::with_capacity(residues * models),
            groups_per_chain: Vec::with_capacity(chains * models),
            chain_ids: Vec::with_capacity(chains * models),
            chain_names: OptionalColumn::Undecided,
        }
    }

    fn append(&mut self, structure: &Structure, group_type_of: &[i32]) -> Result<(), Diagnostic> {
        for chain in structure.data().chains() {
            self.chain_ids
                .push(required_text(chain.label(), "chain label")?);
            self.chain_names.push(
                chain.auth_label().map(str::to_owned),
                "MMTF cannot represent partially absent author chain names",
            )?;
            let residues: Vec<_> = chain.residues().collect();
            self.groups_per_chain
                .push(i32_of(residues.len(), "groups per chain")?);
            for residue in residues {
                self.group_ids.push(residue.auth_seq_id().ok_or_else(|| {
                    refusal("MMTF requires an author residue sequence identifier")
                })?);
                // The dictionary index, not the residue index: identical
                // chemistries share one entry.
                let group_type = group_type_of
                    .get(residue.index().as_usize())
                    .copied()
                    .ok_or_else(|| refusal("MMTF group type is missing for a residue"))?;
                self.group_types.push(group_type);
                self.insertion_codes
                    .push(single_ascii(residue.ins_code(), "insertion code")?);
                self.sequence_indices
                    .push(residue.label_seq_id().map_or(-1, |value| value - 1));
            }
        }
        Ok(())
    }
}

struct AtomColumns {
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    b_factor: OptionalColumn<f32>,
    occupancy: OptionalColumn<f32>,
    atom_id: OptionalColumn<i32>,
    alt_loc: Vec<u8>,
}

impl AtomColumns {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            z: Vec::with_capacity(capacity),
            b_factor: OptionalColumn::Undecided,
            occupancy: OptionalColumn::Undecided,
            atom_id: OptionalColumn::Undecided,
            alt_loc: Vec::with_capacity(capacity),
        }
    }

    fn append(&mut self, structure: &Structure, model: ModelIndex) -> Result<(), Diagnostic> {
        let positions = structure
            .model_positions(model)
            .ok_or_else(|| refusal("MMTF requires a dense coordinate frame"))?;
        if positions.len() != structure.atom_count() as usize {
            return Err(refusal(
                "MMTF frame length differs from the shared topology",
            ));
        }
        for (atom, position) in structure.data().atoms().zip(positions) {
            self.x.push(position[0]);
            self.y.push(position[1]);
            self.z.push(position[2]);
            self.b_factor.push(
                atom.b_factor(),
                "MMTF cannot represent partially absent B factors",
            )?;
            self.occupancy.push(
                atom.occupancy(),
                "MMTF cannot represent partially absent occupancies",
            )?;
            let atom_id = atom
                .atom_site_id()
                .map(|value| {
                    i32::try_from(value)
                        .map_err(|_| refusal("MMTF atom IDs are signed 32-bit integers"))
                })
                .transpose()?;
            self.atom_id.push(
                atom_id,
                "MMTF cannot represent partially absent atom identifiers",
            )?;
            self.alt_loc
                .push(single_ascii(atom.alt_label(), "alternate location")?);
        }
        Ok(())
    }
}

enum OptionalColumn<T> {
    Undecided,
    Present(Vec<T>),
    Absent,
}

impl<T> OptionalColumn<T> {
    fn push(&mut self, value: Option<T>, mismatch: &'static str) -> Result<(), Diagnostic> {
        match (&mut *self, value) {
            (Self::Undecided, Some(value)) => *self = Self::Present(vec![value]),
            (Self::Undecided, None) => *self = Self::Absent,
            (Self::Present(values), Some(value)) => values.push(value),
            (Self::Absent, None) => {}
            (Self::Present(_), None) | (Self::Absent, Some(_)) => return Err(refusal(mismatch)),
        }
        Ok(())
    }

    fn values(&self) -> Option<&[T]> {
        match self {
            Self::Present(values) => Some(values),
            Self::Undecided | Self::Absent => None,
        }
    }
}

type EncodedBonds = (Vec<i32>, Vec<i32>, Vec<i32>);

fn bonds(structure: &Structure, models: usize) -> Result<EncodedBonds, Diagnostic> {
    let count = structure.atom_count() as usize;
    let mut atoms = Vec::with_capacity(structure.data().bonds.len() * models * 2);
    let mut orders = Vec::with_capacity(structure.data().bonds.len() * models);
    let mut resonance = Vec::with_capacity(structure.data().bonds.len() * models);
    for model in 0..models {
        let offset = model
            .checked_mul(count)
            .ok_or_else(|| refusal("MMTF bond offset overflows"))?;
        for bond in structure.data().bonds.iter() {
            let endpoint = |atom: usize| {
                offset
                    .checked_add(atom)
                    .ok_or_else(|| refusal("MMTF bond endpoint overflows"))
                    .and_then(|value| i32_of(value, "bond endpoint"))
            };
            atoms.push(endpoint(bond.atom_a.as_usize())?);
            atoms.push(endpoint(bond.atom_b.as_usize())?);
            let (order, resonant) = encoded_order(bond.order);
            orders.push(order);
            resonance.push(resonant);
        }
    }
    Ok((atoms, orders, resonance))
}

fn entities(
    structure: &Structure,
    metadata: &MmtfMetadata,
    models: usize,
) -> Result<Vec<Entity>, Diagnostic> {
    let mut chains: BTreeMap<EntityIndex, Vec<i32>> = BTreeMap::new();
    let chain_count = structure.chain_count();
    for model in 0..models {
        for chain in structure.data().chains() {
            let entity = chain
                .entity()
                .ok_or_else(|| refusal("chain has no entity"))?;
            let position = model
                .checked_mul(chain_count)
                .and_then(|offset| offset.checked_add(chain.index().as_usize()))
                .ok_or_else(|| refusal("MMTF entity chain index overflows"))?;
            chains
                .entry(entity)
                .or_default()
                .push(i32_of(position, "entity chain index")?);
        }
    }
    structure
        .data()
        .topology
        .entities
        .iter()
        .enumerate()
        .map(|(position, entity)| {
            let source = metadata
                .entities
                .get(position)
                .ok_or_else(|| refusal("MMTF entity metadata is incomplete"))?;
            Ok(Entity {
                chain_index_list: chains
                    .remove(&entity)
                    .ok_or_else(|| refusal("MMTF entity is not instantiated by any chain"))?,
                description: source.description.to_string(),
                kind: source.kind.to_string(),
                sequence: source.sequence.to_string(),
            })
        })
        .collect()
}

include!("writer/helpers.rs");
