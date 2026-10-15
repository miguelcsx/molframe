//! Group-local bond lowering and group lookups for the MMTF reader.
//!
//! Complexity: O(bonds) time per model; O(1) extra space apart from the bonds produced.

use super::Decoded;
use super::bonds::GlobalBond;
use super::validation::{schema_error, usize_of};
use crate::mmtf::schema::Group;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::index::AtomIndex;
use molframe_core::optional::OptionalSymbol;
use molframe_core::{BondOrder, BondProvenance, BondRecord, BondTableBuilder, Element};

/// Adds bonds declared within a single MMTF group.
///
/// Bond endpoints are translated from group-local indices into model-local,
/// retained normalized atom indices. Bonds touching filtered atoms are ignored.
pub(super) fn add_group_bonds(
    group: &Group,
    atom_offset: usize,
    keep_map: &[Option<AtomIndex>],
    bonds: &mut BondTableBuilder,
) -> Result<(), Diagnostic> {
    if !group.bond_atom_list.len().is_multiple_of(2) {
        return Err(schema_error("group bond endpoint list has odd length"));
    }

    for (bond_index, pair) in group.bond_atom_list.as_chunks::<2>().0.iter().enumerate() {
        let first = checked_atom_offset(atom_offset, pair[0])?;
        let second = checked_atom_offset(atom_offset, pair[1])?;

        let Some((atom_a, atom_b)) = retained_atoms(keep_map, first, second) else {
            continue;
        };

        bonds.push(BondRecord {
            atom_a,
            atom_b,
            order: decoded_order(
                &group.bond_order_list,
                &group.bond_resonance_list,
                bond_index,
            ),
            provenance: BondProvenance::File,
            // MMTF carries no `_struct_conn.conn_type_id`; the bond's own order
            // and endpoints are all it says.
            kind: OptionalSymbol::NONE,
        });
    }

    Ok(())
}

/// Adds global MMTF bonds belonging to the selected model.
///
/// Global atom indices outside the model are ignored. Indices inside the model
/// are translated to the model-local retention map, so memory is proportional
/// to one model rather than the complete MMTF file. Bonds touching atoms that
/// were filtered from the normalized structure are also ignored.
pub(super) fn add_global_bonds(
    global_bonds: &[GlobalBond],
    keep_map: &[Option<AtomIndex>],
    bonds: &mut BondTableBuilder,
) -> Result<(), Diagnostic> {
    for bond in global_bonds {
        let first = usize::try_from(bond.atom_a)
            .map_err(|_| schema_error("MMTF local bond endpoint exceeds usize"))?;
        let second = usize::try_from(bond.atom_b)
            .map_err(|_| schema_error("MMTF local bond endpoint exceeds usize"))?;
        let Some((atom_a, atom_b)) = retained_atoms(keep_map, first, second) else {
            continue;
        };

        bonds.push(BondRecord {
            atom_a,
            atom_b,
            order: bond.order,
            provenance: BondProvenance::File,
            // MMTF carries no `_struct_conn.conn_type_id`.
            kind: OptionalSymbol::NONE,
        });
    }

    Ok(())
}

/// Resolves two original atom indices into their retained normalized indices.
///
/// Returns `None` when either endpoint lies outside the mapping or corresponds
/// to an atom that was intentionally filtered out.
pub(super) fn retained_atoms(
    keep_map: &[Option<AtomIndex>],
    first: usize,
    second: usize,
) -> Option<(AtomIndex, AtomIndex)> {
    match (keep_map.get(first), keep_map.get(second)) {
        (Some(Some(atom_a)), Some(Some(atom_b))) => Some((*atom_a, *atom_b)),
        _ => None,
    }
}

/// Converts a group-local signed atom index into an absolute atom index.
///
/// Rejects negative MMTF indices and arithmetic overflow while applying the
/// group atom offset.
pub(super) fn checked_atom_offset(atom_offset: usize, local: i32) -> Result<usize, Diagnostic> {
    atom_offset
        .checked_add(usize_of(local)?)
        .ok_or_else(|| schema_error("MMTF atom index overflows"))
}

/// Decodes MMTF bond-order and resonance fields into the normalized bond order.
///
/// Resonance takes precedence and is represented as aromatic. Missing or
/// unsupported numeric bond orders are represented as unknown.
pub(super) fn decoded_order(orders: &[i32], resonance: &[i32], index: usize) -> BondOrder {
    if resonance.get(index) == Some(&1) {
        return BondOrder::Aromatic;
    }

    match orders.get(index).copied() {
        Some(1) => BondOrder::Single,
        Some(2) => BondOrder::Double,
        Some(3) => BondOrder::Triple,
        Some(4) => BondOrder::Quadruple,
        _ => BondOrder::Unknown,
    }
}

/// Resolves the group template referenced by one decoded group instance.
///
/// Validates both the group-instance index and the referenced group-template
/// index before returning the template.
pub(super) fn group_at(decoded: &Decoded, group_index: usize) -> Result<&Group, Diagnostic> {
    let type_index = decoded
        .group_type
        .get(group_index)
        .copied()
        .ok_or_else(|| schema_error("group type missing"))
        .and_then(usize_of)?;

    decoded
        .file
        .group_list
        .get(type_index)
        .ok_or_else(|| schema_error("group type index is out of bounds"))
}

/// Resolves the chemical element of a group-local atom.
///
/// Missing, out-of-range, or unrecognized element symbols map to
/// `Element::UNKNOWN`.
pub(super) fn group_element(group: &Group, local: usize) -> Element {
    match group
        .element_list
        .as_ref()
        .and_then(|elements| elements.get(local))
        .and_then(|symbol| Element::from_symbol(symbol))
    {
        Some(element) => element,
        None => Element::UNKNOWN,
    }
}

/// Counts atoms contained in a contiguous range of decoded MMTF groups.
///
/// Validates the requested group range, every referenced group template, and
/// all arithmetic used to accumulate the atom count.
pub(super) fn atoms_in_groups(
    decoded: &Decoded,
    start: usize,
    count: usize,
) -> Result<usize, Diagnostic> {
    let end = start
        .checked_add(count)
        .ok_or_else(|| schema_error("MMTF group range overflows"))?;

    let mut atoms = 0usize;

    for group_index in start..end {
        let group_atoms = group_at(decoded, group_index)?.atom_name_list.len();
        atoms = atoms
            .checked_add(group_atoms)
            .ok_or_else(|| schema_error("MMTF atom count overflows"))?;
    }

    Ok(atoms)
}
