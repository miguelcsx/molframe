//! Default bond perception for parsed structures.
//!
//! File connectivity is copied first. Missing edges are then added from
//! polymer attachment rules and a uniform sorted-cell spatial grid. The grid
//! visits 27 neighbouring cells, so the cost is O(N*k) for bounded local
//! density and the output is O(N+B). It retains no per-atom neighbour vectors.

use crate::standard_bonds::standard_bond_order;
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::AtomIndex;
use molframe_core::structure::{AtomRef, Structure};
use num_traits::ToPrimitive;
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct CellKey([i32; 3]);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GridEntry {
    cell: CellKey,
    atom: u32,
}

fn add_grid_bonds(
    structure: &Structure,
    thresholds: &[Option<f32>],
    cutoff: f32,
    output: &mut BondTableBuilder,
    existing: &mut HashSet<(u32, u32)>,
) {
    let positions = structure.positions();
    let mut entries = Vec::with_capacity(positions.len());
    for (atom, position) in positions.iter().enumerate() {
        let Some(cell) = cell_for(*position, cutoff) else {
            continue;
        };
        let Ok(atom) = u32::try_from(atom) else {
            continue;
        };
        entries.push(GridEntry { cell, atom });
    }
    entries.sort_unstable();
    let count = structure.atom_count();
    for first in 0..count {
        let Some(position) = positions.get(first as usize).copied() else {
            continue;
        };
        let Some(base) = cell_for(position, cutoff) else {
            continue;
        };
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let cell = CellKey([
                        base.0[0].saturating_add(dx),
                        base.0[1].saturating_add(dy),
                        base.0[2].saturating_add(dz),
                    ]);
                    let range = cell_range(&entries, cell);
                    for entry in &entries[range] {
                        if entry.atom <= first {
                            continue;
                        }
                        let Some(second_position) = positions.get(entry.atom as usize).copied()
                        else {
                            continue;
                        };
                        let distance_squared = squared_distance(position, second_position);
                        let (Some(left), Some(right)) = (
                            structure.atom(AtomIndex::new(first)),
                            structure.atom(AtomIndex::new(entry.atom)),
                        ) else {
                            continue;
                        };
                        if !candidate_allowed(left, right, thresholds, distance_squared) {
                            continue;
                        }
                        let key = endpoints(first, entry.atom);
                        if !existing.insert(key) {
                            continue;
                        }
                        output.push(BondRecord {
                            atom_a: AtomIndex::new(first),
                            atom_b: AtomIndex::new(entry.atom),
                            order: intra_residue_order(left, right),
                            provenance: BondProvenance::InferredDistance,
                        });
                    }
                }
            }
        }
    }
}

fn cell_range(entries: &[GridEntry], cell: CellKey) -> std::ops::Range<usize> {
    let start = entries.partition_point(|entry| entry.cell < cell);
    let end = entries.partition_point(|entry| entry.cell <= cell);
    start..end
}

fn cell_for(position: [f32; 3], size: f32) -> Option<CellKey> {
    if position.iter().all(|value| value.is_finite()) && size.is_finite() && size > 0.0 {
        Some(CellKey([
            (position[0] / size).floor().to_i32()?,
            (position[1] / size).floor().to_i32()?,
            (position[2] / size).floor().to_i32()?,
        ]))
    } else {
        None
    }
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

fn candidate_allowed(
    left: AtomRef<'_>,
    right: AtomRef<'_>,
    thresholds: &[Option<f32>],
    distance_squared: f32,
) -> bool {
    if left
        .element()
        .is_some_and(molframe_core::Element::is_hydrogen)
        && right
            .element()
            .is_some_and(molframe_core::Element::is_hydrogen)
    {
        return false;
    }
    if !alt_compatible(left, right) {
        return false;
    }
    let (Some(left_threshold), Some(right_threshold)) = (
        thresholds.get(left.index().as_usize()).copied().flatten(),
        thresholds.get(right.index().as_usize()).copied().flatten(),
    ) else {
        return false;
    };
    let maximum = (left_threshold + right_threshold) / 1.95;
    distance_squared >= 0.16 && distance_squared <= maximum * maximum
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
