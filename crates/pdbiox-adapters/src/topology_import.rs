//! Validated native construction from a columnar topology batch.

use crate::{MISSING_STRING, TopologyBatch};
use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
use pdbiox_core::column::Presence;
use pdbiox_core::structure::{CoordinateStore, Structure, StructureData};
use pdbiox_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};
use pdbiox_core::{
    AltId, AtomIndex, BondProvenance, BondRecord, BondTableBuilder, CapacityError, DictionaryFull,
    Element, OptionalI32, OptionalSymbol, ResidueIndex, TableError,
};
use thiserror::Error;

/// Why a columnar topology batch cannot be represented natively.
#[derive(Debug, Error)]
pub enum TopologyImportError {
    /// A column has a different row count from its table.
    #[error("{column} has {observed} rows, expected {expected}")]
    ColumnLength {
        /// Misaligned column.
        column: &'static str,
        /// Rows found.
        observed: usize,
        /// Rows required.
        expected: usize,
    },
    /// An offset column is not monotonic or does not tile its child table.
    #[error("{kind} offset {index} is {observed}, expected at least {minimum}")]
    InvalidOffset {
        /// Parent table.
        kind: &'static str,
        /// Offset position.
        index: usize,
        /// Offset found.
        observed: u32,
        /// Smallest valid offset.
        minimum: u32,
    },
    /// An index exceeds the transferred row count.
    #[error("{kind} index {index} is outside {len} rows")]
    IndexOutOfBounds {
        /// Index domain.
        kind: &'static str,
        /// Invalid index.
        index: usize,
        /// Domain length.
        len: usize,
    },
    /// An atomic number is outside the periodic table.
    #[error("atom {atom} has unsupported atomic number {atomic_number}")]
    InvalidAtomicNumber {
        /// Atom row.
        atom: usize,
        /// Invalid value.
        atomic_number: u8,
    },
    /// The structure-local identifier dictionary is full.
    #[error(transparent)]
    Dictionary(#[from] DictionaryFull),
    /// A compact hierarchy table rejected a row.
    #[error(transparent)]
    Table(#[from] TableError),
    /// An entity table exceeded its compact index space.
    #[error(transparent)]
    EntityCapacity(#[from] CapacityError),
}

impl TopologyBatch {
    /// Consumes the columns and builds one immutable native model.
    ///
    /// # Errors
    ///
    /// Returns an explicit error when columns are misaligned, hierarchy
    /// offsets or indices are invalid, or native compact storage is exhausted.
    pub fn into_structure(self) -> Result<Structure, TopologyImportError> {
        validate_columns(&self)?;
        build_structure(&self)
    }
}

fn validate_columns(batch: &TopologyBatch) -> Result<(), TopologyImportError> {
    let chains = batch.chain_ids.len();
    let residues = batch.residue_names.len();
    let atoms = batch.atom_names.len();
    let bonds = batch.bond_atom_a.len();
    check_len(
        "chain_residue_offsets",
        batch.chain_residue_offsets.len(),
        chains + 1,
    )?;
    for (name, observed) in [
        ("residue_numbers", batch.residue_numbers.len()),
        (
            "residue_insertion_codes",
            batch.residue_insertion_codes.len(),
        ),
        ("residue_chain", batch.residue_chain.len()),
    ] {
        check_len(name, observed, residues)?;
    }
    check_len(
        "residue_number_validity",
        batch.residue_number_validity.len() as usize,
        residues,
    )?;
    check_len(
        "residue_is_heterogen",
        batch.residue_is_heterogen.len() as usize,
        residues,
    )?;
    check_len(
        "residue_atom_offsets",
        batch.residue_atom_offsets.len(),
        residues + 1,
    )?;
    for (name, observed) in [
        ("atomic_numbers", batch.atomic_numbers.len()),
        ("masses", batch.masses.len()),
        ("atom_serials", batch.atom_serials.len()),
        ("formal_charges", batch.formal_charges.len()),
        ("occupancies", batch.occupancies.len()),
        ("b_factors", batch.b_factors.len()),
        (
            "atom_alternate_locations",
            batch.atom_alternate_locations.len(),
        ),
        ("atom_residue", batch.atom_residue.len()),
        ("position_x", batch.position_x.len()),
        ("position_y", batch.position_y.len()),
        ("position_z", batch.position_z.len()),
    ] {
        check_len(name, observed, atoms)?;
    }
    for (name, observed) in [
        ("atom_serial_validity", batch.atom_serial_validity.len()),
        ("formal_charge_validity", batch.formal_charge_validity.len()),
        ("occupancy_validity", batch.occupancy_validity.len()),
        ("b_factor_validity", batch.b_factor_validity.len()),
    ] {
        check_len(name, observed as usize, atoms)?;
    }
    check_len("bond_atom_b", batch.bond_atom_b.len(), bonds)?;
    check_len("bond_orders", batch.bond_orders.len(), bonds)?;
    validate_offsets("chain", &batch.chain_residue_offsets, residues)?;
    validate_offsets("residue", &batch.residue_atom_offsets, atoms)?;
    validate_parents(batch)?;
    validate_string_ids(batch)?;
    for (atom, atomic_number) in batch.atomic_numbers.iter().copied().enumerate() {
        if !(1..=118).contains(&atomic_number) {
            return Err(TopologyImportError::InvalidAtomicNumber {
                atom,
                atomic_number,
            });
        }
    }
    for endpoint in batch.bond_atom_a.iter().chain(&batch.bond_atom_b) {
        if *endpoint as usize >= atoms {
            return Err(TopologyImportError::IndexOutOfBounds {
                kind: "bond endpoint",
                index: *endpoint as usize,
                len: atoms,
            });
        }
    }
    Ok(())
}

fn validate_parents(batch: &TopologyBatch) -> Result<(), TopologyImportError> {
    for chain in 0..batch.chain_ids.len() {
        for residue in batch.chain_residue_offsets[chain]..batch.chain_residue_offsets[chain + 1] {
            if batch.residue_chain[residue as usize] as usize != chain {
                return Err(TopologyImportError::IndexOutOfBounds {
                    kind: "residue parent",
                    index: batch.residue_chain[residue as usize] as usize,
                    len: batch.chain_ids.len(),
                });
            }
        }
    }
    for residue in 0..batch.residue_names.len() {
        for atom in batch.residue_atom_offsets[residue]..batch.residue_atom_offsets[residue + 1] {
            if batch.atom_residue[atom as usize] as usize != residue {
                return Err(TopologyImportError::IndexOutOfBounds {
                    kind: "atom parent",
                    index: batch.atom_residue[atom as usize] as usize,
                    len: batch.residue_names.len(),
                });
            }
        }
    }
    Ok(())
}

fn validate_string_ids(batch: &TopologyBatch) -> Result<(), TopologyImportError> {
    for id in batch
        .chain_ids
        .iter()
        .chain(&batch.residue_names)
        .chain(&batch.atom_names)
    {
        string_at(batch, *id, "required string")?;
    }
    for id in batch
        .residue_insertion_codes
        .iter()
        .chain(&batch.atom_alternate_locations)
        .filter(|id| **id != MISSING_STRING)
    {
        string_at(batch, *id, "optional string")?;
    }
    Ok(())
}

fn check_len(
    column: &'static str,
    observed: usize,
    expected: usize,
) -> Result<(), TopologyImportError> {
    if observed == expected {
        Ok(())
    } else {
        Err(TopologyImportError::ColumnLength {
            column,
            observed,
            expected,
        })
    }
}

fn validate_offsets(
    kind: &'static str,
    offsets: &[u32],
    child_count: usize,
) -> Result<(), TopologyImportError> {
    let mut previous = 0;
    for (index, observed) in offsets.iter().copied().enumerate() {
        if observed < previous {
            return Err(TopologyImportError::InvalidOffset {
                kind,
                index,
                observed,
                minimum: previous,
            });
        }
        previous = observed;
    }
    if previous as usize != child_count {
        return Err(TopologyImportError::ColumnLength {
            column: kind,
            observed: previous as usize,
            expected: child_count,
        });
    }
    Ok(())
}

fn build_structure(batch: &TopologyBatch) -> Result<Structure, TopologyImportError> {
    let mut data = StructureData::empty();
    let mut builder = ChunkBuilder::new();
    builder.reserve(batch.atom_names.len());
    builder.start_model(0);
    for chain_index in 0..batch.chain_ids.len() {
        let chain_text = string_at(batch, batch.chain_ids[chain_index], "chain string")?;
        let chain_id = data.dictionary.intern(chain_text)?;
        let entity = data.topology.entities.push(
            chain_id,
            EntityKind::Unknown,
            OptionalSymbol::NONE,
            &[],
        )?;
        let residue_start = batch.chain_residue_offsets[chain_index];
        let residue_end = batch.chain_residue_offsets[chain_index + 1];
        for residue in residue_start..residue_end {
            push_residue_and_atoms(&mut data, &mut builder, residue, batch)?;
        }
        data.topology.chains.push(
            ChainRecord {
                label_asym_id: chain_id,
                auth_asym_id: OptionalSymbol::some(chain_id),
                entity,
                polymer_kind: PolymerKind::None,
            },
            residue_start..residue_end,
        )?;
    }
    data.topology.models.push(
        1,
        0..u32::try_from(batch.chain_ids.len()).map_err(|_| TopologyImportError::ColumnLength {
            column: "chains",
            observed: batch.chain_ids.len(),
            expected: u32::MAX as usize,
        })?,
    )?;
    let mut bonds = BondTableBuilder::new();
    for ((atom_a, atom_b), order) in batch
        .bond_atom_a
        .iter()
        .zip(&batch.bond_atom_b)
        .zip(&batch.bond_orders)
    {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(*atom_a),
            atom_b: AtomIndex::new(*atom_b),
            order: *order,
            provenance: BondProvenance::File,
        });
    }
    let (chunks, coordinates) = builder.finish();
    data.chunks = chunks.into();
    data.coords = CoordinateStore::Single(coordinates);
    data.bonds = bonds.finish();
    Ok(Structure::new(data))
}

fn push_residue_and_atoms(
    data: &mut StructureData,
    builder: &mut ChunkBuilder,
    residue: u32,
    batch: &TopologyBatch,
) -> Result<(), TopologyImportError> {
    let index = residue as usize;
    let residue_name = data.dictionary.intern(string_at(
        batch,
        batch.residue_names[index],
        "residue string",
    )?)?;
    let insertion = optional_symbol(
        &mut data.dictionary,
        optional_string(batch, batch.residue_insertion_codes[index])?,
    )?;
    let atom_start = batch.residue_atom_offsets[index];
    let atom_end = batch.residue_atom_offsets[index + 1];
    data.topology.residues.push(
        ResidueRecord {
            label_comp_id: residue_name,
            auth_comp_id: OptionalSymbol::some(residue_name),
            label_seq_id: OptionalI32::from(valid_value(
                &batch.residue_numbers,
                &batch.residue_number_validity,
                residue,
            )),
            auth_seq_id: OptionalI32::from(valid_value(
                &batch.residue_numbers,
                &batch.residue_number_validity,
                residue,
            )),
            ins_code: insertion,
            het: batch.residue_is_heterogen.test(residue),
        },
        atom_start..atom_end,
    )?;
    for atom in atom_start..atom_end {
        push_atom(data, builder, atom, residue, batch)?;
    }
    Ok(())
}

fn push_atom(
    data: &mut StructureData,
    builder: &mut ChunkBuilder,
    atom: u32,
    residue: u32,
    batch: &TopologyBatch,
) -> Result<(), TopologyImportError> {
    let index = atom as usize;
    let atom_name =
        data.dictionary
            .intern(string_at(batch, batch.atom_names[index], "atom string")?)?;
    let alternate = optional_symbol(
        &mut data.dictionary,
        optional_string(batch, batch.atom_alternate_locations[index])?,
    )?;
    let alt_id = match alternate.get().and_then(AltId::labelled) {
        Some(value) => value,
        None => AltId::BLANK,
    };
    let atom_site_id = match valid_value(&batch.atom_serials, &batch.atom_serial_validity, atom) {
        Some(value) => value,
        None => atom
            .checked_add(1)
            .ok_or(TopologyImportError::ColumnLength {
                column: "atom serial",
                observed: atom as usize,
                expected: u32::MAX as usize,
            })?,
    };
    builder.push(AtomRecord {
        position: Some([
            batch.position_x[index],
            batch.position_y[index],
            batch.position_z[index],
        ]),
        element: Element::from_atomic_number(batch.atomic_numbers[index]),
        atom_name,
        auth_atom_name: OptionalSymbol::some(atom_name),
        alternate_component_id: OptionalSymbol::NONE,
        alt_id,
        residue: ResidueIndex::new(residue),
        occupancy: present_value(&batch.occupancies, &batch.occupancy_validity, atom),
        b_factor: present_value(&batch.b_factors, &batch.b_factor_validity, atom),
        formal_charge: present_value(&batch.formal_charges, &batch.formal_charge_validity, atom),
        atom_site_id,
    });
    Ok(())
}

fn string_at<'a>(
    batch: &'a TopologyBatch,
    id: u32,
    kind: &'static str,
) -> Result<&'a str, TopologyImportError> {
    batch
        .string(id)
        .ok_or(TopologyImportError::IndexOutOfBounds {
            kind,
            index: id as usize,
            len: batch.strings.len(),
        })
}

fn optional_string(batch: &TopologyBatch, id: u32) -> Result<Option<&str>, TopologyImportError> {
    if id == MISSING_STRING {
        Ok(None)
    } else {
        string_at(batch, id, "optional string").map(Some)
    }
}

fn optional_symbol(
    dictionary: &mut pdbiox_core::Interner,
    value: Option<&str>,
) -> Result<OptionalSymbol, DictionaryFull> {
    match value {
        Some(text) => dictionary.intern(text).map(OptionalSymbol::some),
        None => Ok(OptionalSymbol::NONE),
    }
}

fn valid_value<T: Copy>(
    values: &[T],
    validity: &pdbiox_core::column::ValidityMask,
    position: u32,
) -> Option<T> {
    if validity.get(position).is_present() {
        values.get(position as usize).copied()
    } else {
        None
    }
}

fn present_value<T: Copy + Default>(
    values: &[T],
    validity: &pdbiox_core::column::ValidityMask,
    position: u32,
) -> (T, Presence) {
    let value = match values.get(position as usize).copied() {
        Some(value) => value,
        None => T::default(),
    };
    (value, validity.get(position))
}
