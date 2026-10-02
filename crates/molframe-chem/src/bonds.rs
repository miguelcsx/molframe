//! Default bond perception for parsed structures.
//!
//! File connectivity is copied first. Missing edges are then added from
//! polymer attachment rules and a uniform sorted-cell spatial grid. The grid
//! visits 27 neighbouring cells, so the cost is O(N*k) for bounded local
//! density and the output is O(N+B). It retains no per-atom neighbour vectors.

use crate::grid::{CellGrid, cell_for};
use crate::standard_bonds::standard_bond_order;
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::AtomIndex;
use molframe_core::structure::{AtomRef, Structure};
use std::collections::HashSet;

/// Adds chemically plausible missing bonds to one immutable structure.
///
/// Existing file records always win. The returned snapshot shares all source
/// storage with the input except for its compact bond columns.
///
/// # Errors
///
/// Returns a diagnostic when coordinates are not a dense first-model column.
pub fn perceive_bonds(structure: &Structure) -> Result<Structure, Diagnostic> {
    if !structure.data().coords.is_dense()
        || structure.positions().len() != structure.atom_count() as usize
    {
        return Err(Diagnostic::new(Code::E6008)
            .with_message("default bond perception requires one dense coordinate column"));
    }
    let mut output = BondTableBuilder::new();
    let mut existing = HashSet::with_capacity(structure.data().bonds.len().saturating_mul(2));
    for bond in structure.data().bonds.iter() {
        existing.insert(endpoints(bond.atom_a.get(), bond.atom_b.get()));
        output.push(bond);
    }

    add_polymer_links(structure, &mut output, &mut existing);
    let thresholds: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| atom.element().and_then(element_threshold))
        .collect();
    let maximum = thresholds.iter().flatten().copied().fold(0.0_f32, f32::max) * 2.0 / 1.95;
    if maximum > 0.0 && maximum.is_finite() {
        add_grid_bonds(structure, &thresholds, maximum, &mut output, &mut existing);
    }

    let mut data = structure.data().clone();
    data.bonds = output.finish();
    Ok(Structure::new(data))
}

/// One atom that can bond, stored where the grid keeps it.
#[derive(Clone, Copy)]
struct GridAtom {
    atom: u32,
    position: [f32; 3],
    threshold: f32,
}

fn add_grid_bonds(
    structure: &Structure,
    thresholds: &[Option<f32>],
    cutoff: f32,
    output: &mut BondTableBuilder,
    existing: &mut HashSet<(u32, u32)>,
) {
    let mut entries = Vec::with_capacity(thresholds.len());
    for (atom, (position, threshold)) in structure.positions().iter().zip(thresholds).enumerate() {
        let (Some(threshold), Some(cell), Ok(atom)) =
            (*threshold, cell_for(*position, cutoff), u32::try_from(atom))
        else {
            continue;
        };
        entries.push((
            cell,
            GridAtom {
                atom,
                position: *position,
                threshold,
            },
        ));
    }
    let grid = CellGrid::build(entries);
    grid.for_each_cell(|own, neighbourhood| {
        for first in &grid.items()[own] {
            for range in neighbourhood {
                for second in &grid.items()[range.clone()] {
                    if second.atom <= first.atom {
                        continue;
                    }
                    // The distance window rejects nearly every candidate, so it
                    // runs before anything that has to look an atom up.
                    let distance_squared = squared_distance(first.position, second.position);
                    if !in_bonding_window(first.threshold, second.threshold, distance_squared) {
                        continue;
                    }
                    let (Some(left), Some(right)) = (
                        structure.atom(AtomIndex::new(first.atom)),
                        structure.atom(AtomIndex::new(second.atom)),
                    ) else {
                        continue;
                    };
                    if !pair_compatible(left, right)
                        || !existing.insert(endpoints(first.atom, second.atom))
                    {
                        continue;
                    }
                    output.push(BondRecord {
                        atom_a: AtomIndex::new(first.atom),
                        atom_b: AtomIndex::new(second.atom),
                        order: intra_residue_order(left, right),
                        provenance: BondProvenance::InferredDistance,
                    });
                }
            }
        }
    });
}

fn add_polymer_links(
    structure: &Structure,
    output: &mut BondTableBuilder,
    existing: &mut HashSet<(u32, u32)>,
) {
    for chain in structure.data().chains() {
        let residues: Vec<_> = chain.residues().collect();
        for pair in residues.windows(2) {
            let [left_residue, right_residue] = pair else {
                continue;
            };
            if let Some((left, right)) = left_residue.atom("C").zip(right_residue.atom("N"))
                && within(left, right, 1.75)
            {
                add_link(left, right, output, existing);
            }
            let left = ["O3'", "O3*", "O3"]
                .iter()
                .find_map(|name| left_residue.atom(name));
            if let Some((left, right)) = left.zip(right_residue.atom("P"))
                && within(left, right, 1.9)
            {
                add_link(left, right, output, existing);
            }
        }
    }
}

fn add_link(
    left: AtomRef<'_>,
    right: AtomRef<'_>,
    output: &mut BondTableBuilder,
    existing: &mut HashSet<(u32, u32)>,
) {
    let key = endpoints(left.index().get(), right.index().get());
    if existing.insert(key) {
        output.push(BondRecord {
            atom_a: left.index(),
            atom_b: right.index(),
            order: BondOrder::Polymeric,
            provenance: BondProvenance::InferredDistance,
        });
    }
}

/// Whether a distance lies between the closest plausible contact and the sum of
/// the two elements' bonding thresholds.
fn in_bonding_window(left_threshold: f32, right_threshold: f32, distance_squared: f32) -> bool {
    let maximum = (left_threshold + right_threshold) / 1.95;
    distance_squared >= 0.16 && distance_squared <= maximum * maximum
}

/// Whether two atoms may bond at all: two hydrogens never do, and two atoms in
/// different alternate locations do not coexist.
fn pair_compatible(left: AtomRef<'_>, right: AtomRef<'_>) -> bool {
    let both_hydrogen = left
        .element()
        .is_some_and(molframe_core::Element::is_hydrogen)
        && right
            .element()
            .is_some_and(molframe_core::Element::is_hydrogen);
    !both_hydrogen && alt_compatible(left, right)
}

fn intra_residue_order(left: AtomRef<'_>, right: AtomRef<'_>) -> BondOrder {
    let (Some(left_residue), Some(right_residue)) = (left.residue(), right.residue()) else {
        return BondOrder::Single;
    };
    if left_residue.index() != right_residue.index() {
        return BondOrder::Single;
    }
    let (Some(component), Some(left_name), Some(right_name)) =
        (left.component_name(), left.name(), right.name())
    else {
        return BondOrder::Single;
    };
    standard_bond_order(component, left_name, right_name)
}

fn alt_compatible(left: AtomRef<'_>, right: AtomRef<'_>) -> bool {
    match (left.alt_id(), right.alt_id()) {
        (Some(left), Some(right)) => left.is_blank() || right.is_blank() || left == right,
        _ => true,
    }
}

fn within(left: AtomRef<'_>, right: AtomRef<'_>, maximum: f32) -> bool {
    let (Some(left), Some(right)) = (left.position(), right.position()) else {
        return false;
    };
    squared_distance(left, right) <= maximum * maximum
}

fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}

fn endpoints(left: u32, right: u32) -> (u32, u32) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

fn element_threshold(element: molframe_core::Element) -> Option<f32> {
    match element.symbol() {
        "H" => Some(1.42),
        "C" => Some(1.75),
        "N" => Some(1.60),
        "O" => Some(1.52),
        "P" => Some(2.00),
        "S" => Some(1.90),
        "F" | "Cl" | "Br" | "I" => Some(1.80),
        symbol if is_metal(symbol) => Some(2.70),
        _ => None,
    }
}

fn is_metal(symbol: &str) -> bool {
    matches!(
        symbol,
        "Li" | "Na"
            | "K"
            | "Mg"
            | "Ca"
            | "Mn"
            | "Fe"
            | "Co"
            | "Ni"
            | "Cu"
            | "Zn"
            | "Mo"
            | "Ag"
            | "Cd"
            | "Hg"
            | "Au"
            | "Al"
            | "Si"
            | "Cr"
            | "V"
            | "W"
    )
}

#[cfg(test)]
#[path = "bonds_tests.rs"]
mod tests;
