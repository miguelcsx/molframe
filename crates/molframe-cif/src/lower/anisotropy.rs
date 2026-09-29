//! Anisotropic displacement lowering from `atom_site_anisotrop`.
//!
//! The ellipsoids are their own mmCIF category, kept beside `atom_site` and
//! joined to atoms by the deposited `_atom_site.id` value. Row order carries no
//! guarantee of agreeing with the coordinate rows, and files that omit the
//! category simply carry no table, so absence leaves `anisotropy` unresolved
//! rather than falsely reported as an empty known set.
//!
//! The legacy inline spelling keeps `aniso_*` items inside `atom_site` itself,
//! and the older `U_11`-style item names are accepted beside the bracketed
//! modern names. The B-tensor form converts on the way in with
//! `U = B / (8π²)`; a row that mixes the two forms inconsistently, or that is
//! missing a component, is left unattached.
//!
//! Cost: one pass over the atoms builds the identifier join and one pass over
//! the anisotropy rows attaches each row, so the total is O(atoms + aniso
//! rows) with a hash map holding only the joined atoms.

use crate::document::{Category, CifValue, DataBlock};
use molframe_core::anisotropy::{AnisotropicDisplacement, AnisotropyTableBuilder};
use molframe_core::diagnostic::Diagnostics;
use molframe_core::index::AtomIndex;
use molframe_core::structure::StructureData;
use num_traits::ToPrimitive;
use std::collections::HashMap;

const CATEGORY: &str = "atom_site_anisotrop";
/// The spelling that keeps the ellipsoid inside `atom_site` itself.
const INLINE: &str = "atom_site";
/// B factors are displacement parameters in ångström squared scaled by 8π², so
/// one B value converts to its U counterpart by dividing by that factor.
const B_TO_U: f64 = 0.012_665_147_955_292_222;

/// One tensor component, with every item spelling the file may have used for
/// it under either form.
struct Component {
    u: &'static [&'static str],
    b: &'static [&'static str],
}

const U_11: &[&str] = &["U[1][1]", "u[1][1]", "U_11", "u_11"];
const U_22: &[&str] = &["U[2][2]", "u[2][2]", "U_22", "u_22"];
const U_33: &[&str] = &["U[3][3]", "u[3][3]", "U_33", "u_33"];
const U_12: &[&str] = &["U[1][2]", "u[1][2]", "U_12", "u_12"];
const U_13: &[&str] = &["U[1][3]", "u[1][3]", "U_13", "u_13"];
const U_23: &[&str] = &["U[2][3]", "u[2][3]", "U_23", "u_23"];
const B_11: &[&str] = &["B[1][1]", "b[1][1]", "B_11", "b_11"];
const B_22: &[&str] = &["B[2][2]", "b[2][2]", "B_22", "b_22"];
const B_33: &[&str] = &["B[3][3]", "b[3][3]", "B_33", "b_33"];
const B_12: &[&str] = &["B[1][2]", "b[1][2]", "B_12", "b_12"];
const B_13: &[&str] = &["B[1][3]", "b[1][3]", "B_13", "b_13"];
const B_23: &[&str] = &["B[2][3]", "b[2][3]", "B_23", "b_23"];

/// The dedicated category's spellings, in `[U11, U22, U33, U12, U13, U23]`
/// order.
const CATEGORY_COMPONENTS: [Component; 6] = [
    Component { u: U_11, b: B_11 },
    Component { u: U_22, b: B_22 },
    Component { u: U_33, b: B_33 },
    Component { u: U_12, b: B_12 },
    Component { u: U_13, b: B_13 },
    Component { u: U_23, b: B_23 },
];
/// The inline spelling prefixes the same item names with `aniso_`.
const INLINE_COMPONENTS: [Component; 6] = [
    Component {
        u: &["aniso_U[1][1]", "aniso_u[1][1]"],
        b: &["aniso_B[1][1]", "aniso_b[1][1]"],
    },
    Component {
        u: &["aniso_U[2][2]", "aniso_u[2][2]"],
        b: &["aniso_B[2][2]", "aniso_b[2][2]"],
    },
    Component {
        u: &["aniso_U[3][3]", "aniso_u[3][3]"],
        b: &["aniso_B[3][3]", "aniso_b[3][3]"],
    },
    Component {
        u: &["aniso_U[1][2]", "aniso_u[1][2]"],
        b: &["aniso_B[1][2]", "aniso_b[1][2]"],
    },
    Component {
        u: &["aniso_U[1][3]", "aniso_u[1][3]"],
        b: &["aniso_B[1][3]", "aniso_b[1][3]"],
    },
    Component {
        u: &["aniso_U[2][3]", "aniso_u[2][3]"],
        b: &["aniso_B[2][3]", "aniso_b[2][3]"],
    },
];

pub(super) fn read(block: &DataBlock, data: &mut StructureData, _findings: &mut Diagnostics) {
    let Some((category, components)) = select_source(block) else {
        return;
    };
    let atoms_by_id = atoms_by_site_id(data);
    if atoms_by_id.is_empty() {
        // No deposited identifiers to join on: attaching by row position would
        // be a guess, so the table stays absent.
        return;
    }
    let mut builder = AnisotropyTableBuilder::new();
    for row in 0..category.row_count() {
        let Some(atom) = source_atom(category, row, &atoms_by_id) else {
            continue;
        };
        let Some(tensor) = source_tensor(category, row, components) else {
            continue;
        };
        builder.push(AnisotropicDisplacement { atom, u: tensor });
    }
    data.anisotropy = builder.finish();
}

/// Picks the first source the block actually carries, dedicated category first.
fn select_source(block: &DataBlock) -> Option<(&Category, &'static [Component; 6])> {
    for (category, components) in [
        (block.category(CATEGORY), &CATEGORY_COMPONENTS),
        (block.category(INLINE), &INLINE_COMPONENTS),
    ] {
        if let Some(category) = category {
            return Some((category, components));
        }
    }
    None
}

/// Maps every atom by its deposited `_atom_site.id`.
fn atoms_by_site_id(data: &StructureData) -> HashMap<u32, AtomIndex> {
    let mut map = HashMap::new();
    for chunk in data.chunks.iter() {
        let range = chunk.atoms();
        for atom in range.start..range.end {
            let local = atom - range.start;
            if let Some(id) = chunk.atom_site_id(local) {
                map.entry(id).or_insert(AtomIndex::new(atom));
            }
        }
    }
    map
}

fn source_atom(
    category: &Category,
    row: usize,
    atoms_by_id: &HashMap<u32, AtomIndex>,
) -> Option<AtomIndex> {
    let joined = category.value("id", row)?.as_integer()?;
    let id = u32::try_from(joined).ok()?;
    atoms_by_id.get(&id).copied()
}

/// Reads one row's tensor, preferring the U form and converting the B form.
///
/// A component whose every spelling is absent, or present as a sentinel, makes
/// the row unattachable: guessing a diagonal-only tensor would invent data the
/// file never stated.
fn source_tensor(category: &Category, row: usize, components: &[Component; 6]) -> Option<[f32; 6]> {
    let mut u = [0.0_f32; 6];
    let mut b = [0.0_f32; 6];
    let mut u_present = 0;
    let mut b_present = 0;
    for (position, component) in components.iter().enumerate() {
        if let Some(value) = component_value(category, row, component.u) {
            u[position] = value;
            u_present += 1;
        }
        if let Some(value) = component_value(category, row, component.b) {
            b[position] = value;
            b_present += 1;
        }
    }
    if u_present == 6 {
        return Some(u);
    }
    if b_present == 6 {
        return b
            .iter()
            .map(|value| b_to_u(*value))
            .collect::<Option<Vec<_>>>()
            .map(|values| {
                [
                    values[0], values[1], values[2], values[3], values[4], values[5],
                ]
            });
    }
    None
}

/// Converts one B-form component to its U counterpart, ÷ 8π².
fn b_to_u(value: f32) -> Option<f32> {
    (f64::from(value) * B_TO_U).to_f32()
}

fn component_value(category: &Category, row: usize, spellings: &[&str]) -> Option<f32> {
    for item in spellings {
        if let Some(value) = category.value(item, row).and_then(CifValue::as_float) {
            return value.to_f32();
        }
    }
    None
}

#[cfg(test)]
#[path = "anisotropy_tests.rs"]
mod tests;
