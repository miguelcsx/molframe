//! Atom, residue and entity record construction for the MMTF reader.
//!
//! Complexity: O(1) per record; O(chains) for the entity map.

use super::Decoded;
use super::group_bonds::group_element;
use super::validation::{
    char_symbol, entity_kind, intern, optional_float, position, schema_error, usize_of,
};
use crate::mmtf::schema::Group;
use molframe_core::chunk::AtomRecord;
use molframe_core::column::Presence;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::index::{EntityIndex, ResidueIndex};
use molframe_core::optional::{OptionalI32, OptionalSymbol};
use molframe_core::structure::StructureData;
use molframe_core::symbol::AltId;
use molframe_core::topology::{EntityKind, ResidueRecord};

/// Builds one normalized atom record from an MMTF group entry.
///
/// Resolves interned identifiers, coordinates, alternate locations, optional
/// experimental values, formal charge, and the source atom-site identifier.
pub(super) fn atom_record(
    data: &mut StructureData,
    decoded: &Decoded,
    group: &Group,
    local: usize,
    atom: usize,
    residue: ResidueIndex,
) -> Result<AtomRecord, Diagnostic> {
    let atom_name = intern(
        data,
        group
            .atom_name_list
            .get(local)
            .ok_or_else(|| schema_error("group atom name missing"))?,
    )?;

    let alt_symbol = char_symbol(
        data,
        decoded
            .alt_loc
            .as_ref()
            .and_then(|values| values.get(atom))
            .copied(),
    )?;

    let alt_id = match alt_symbol {
        Some(symbol) => AltId::labelled(symbol)
            .ok_or_else(|| schema_error("alternate-location identifier exceeds its encoding"))?,
        None => AltId::BLANK,
    };

    let charge = group
        .formal_charge_list
        .get(local)
        .copied()
        .ok_or_else(|| schema_error("group formal charge missing"))?;

    let formal_charge =
        i8::try_from(charge).map_err(|_| schema_error("formal charge is outside i8"))?;

    Ok(AtomRecord {
        position: Some(position(decoded, atom)?),
        element: group_element(group, local),
        atom_name,
        auth_atom_name: OptionalSymbol::NONE,
        alternate_component_id: OptionalSymbol::NONE,
        alt_id,
        residue,
        occupancy: optional_float(decoded.occupancy.as_ref(), atom, 1.0),
        b_factor: optional_float(decoded.b_factor.as_ref(), atom, 0.0),
        formal_charge: (formal_charge, Presence::Present),
        atom_site_id: atom_site_id(decoded.atom_id.as_ref(), atom)?,
    })
}

/// Builds one normalized residue record from an MMTF group.
///
/// Converts component identifiers, sequence identifiers, insertion codes, and
/// entity classification into the normalized hierarchy representation.
pub(super) fn residue_record(
    data: &mut StructureData,
    decoded: &Decoded,
    group: &Group,
    group_index: usize,
    kind: EntityKind,
) -> Result<ResidueRecord, Diagnostic> {
    let component = intern(data, &group.name)?;

    let sequence = decoded
        .sequence_index
        .as_ref()
        .and_then(|values| values.get(group_index))
        .and_then(|value| value.checked_add(1))
        .map_or(OptionalI32::NONE, OptionalI32::some);

    let auth = decoded
        .group_id
        .get(group_index)
        .copied()
        .ok_or_else(|| schema_error("group id missing"))?;

    let ins_code = char_symbol(
        data,
        decoded
            .ins_code
            .as_ref()
            .and_then(|values| values.get(group_index))
            .copied(),
    )?
    .map_or(OptionalSymbol::NONE, OptionalSymbol::some);

    Ok(ResidueRecord {
        label_comp_id: component,
        auth_comp_id: OptionalSymbol::some(component),
        label_seq_id: sequence,
        auth_seq_id: OptionalI32::some(auth),
        ins_code,
        het: kind != EntityKind::Polymer,
    })
}

/// Builds normalized entities and maps every MMTF chain to its owning entity.
///
/// Rejects missing entity metadata, invalid chain indices, duplicate ownership,
/// and chains that are not assigned to any entity.
pub(super) fn build_entities(
    data: &mut StructureData,
    decoded: &Decoded,
) -> Result<Vec<EntityIndex>, Diagnostic> {
    let chain_count = usize_of(decoded.file.num_chains)?;
    let entities =
        decoded.file.entity_list.as_deref().ok_or_else(|| {
            schema_error("MMTF entityList is required by the normalized hierarchy")
        })?;

    let mut mapping = vec![None; chain_count];

    for (position, entity) in entities.iter().enumerate() {
        let ordinal = position
            .checked_add(1)
            .ok_or_else(|| schema_error("MMTF entity ordinal overflows"))?;
        let ordinal = ordinal.to_string();
        let id = intern(data, &ordinal)?;

        let description = if entity.description.is_empty() {
            OptionalSymbol::NONE
        } else {
            OptionalSymbol::some(intern(data, &entity.description)?)
        };

        let index = data
            .topology
            .entities
            .push(id, entity_kind(entity)?, description, &[])
            .map_err(|error| {
                schema_error("entity table rejected a row").with_context("cause", error.to_string())
            })?;

        for &chain in &entity.chain_index_list {
            let chain = usize_of(chain)?;
            let slot = mapping
                .get_mut(chain)
                .ok_or_else(|| schema_error("entity chain index is out of bounds"))?;

            if slot.replace(index).is_some() {
                return Err(schema_error("chain belongs to more than one MMTF entity"));
            }
        }
    }

    mapping
        .into_iter()
        .map(|entity| entity.ok_or_else(|| schema_error("entity mapping missing")))
        .collect()
}

/// Reads and validates an optional MMTF atom-site identifier.
///
/// An absent identifier list is represented by zero. When the source list is
/// present, identifiers must exist, be non-negative, and be non-zero.
pub(super) fn atom_site_id(values: Option<&Vec<i32>>, atom: usize) -> Result<u32, Diagnostic> {
    let Some(values) = values else {
        return Ok(0);
    };

    let value = values
        .get(atom)
        .copied()
        .ok_or_else(|| schema_error("atom identifier missing"))?;

    let value = u32::try_from(value).map_err(|_| schema_error("atom identifier is negative"))?;

    if value == 0 {
        Err(schema_error(
            "present MMTF atom identifier cannot be represented as zero",
        ))
    } else {
        Ok(value)
    }
}
