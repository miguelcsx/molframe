//! Default bond perception for parsed structures.
//!
//! File connectivity is copied first. Missing edges are then added from
//! polymer attachment rules and a uniform sorted-cell spatial grid. The grid
//! visits 27 neighbouring cells, so the cost is O(N*k) for bounded local
//! density and the output is O(N+B). It retains no per-atom neighbour vectors.

use crate::grid::{CellGrid, cell_for};
use crate::standard_bonds::standard_bond_order;
use molframe_core::ExecutionContext;
use molframe_core::bond::{BondOrder, BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::AtomIndex;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use molframe_core::structure::{AtomRef, Structure};
use std::ops::Range;

/// Adds chemically plausible missing bonds to one immutable structure.
///
/// Existing file records always win. The returned snapshot shares all source
/// storage with the input except for its compact bond columns.
///
/// # Errors
///
/// Returns a diagnostic when coordinates are not a dense first-model column.
pub fn perceive_bonds(structure: &Structure) -> Result<Structure, Diagnostic> {
    perceive_bonds_in(structure, &ExecutionContext::default())
}

/// [`perceive_bonds`] on the worker budget of `context`.
///
/// The result is identical at every worker count: the search is divided into
/// fixed cell blocks whose bonds are merged in block order.
///
/// # Errors
///
/// Returns a diagnostic when coordinates are not a dense first-model column,
/// or when a worker thread panics.
pub fn perceive_bonds_in(
    structure: &Structure,
    context: &ExecutionContext,
) -> Result<Structure, Diagnostic> {
    if !structure.data().coords.is_dense()
        || structure.positions().len() != structure.atom_count() as usize
    {
        return Err(Diagnostic::new(Code::E6008)
            .with_message("default bond perception requires one dense coordinate column"));
    }
    // Provenance is decided by insertion order: `BondTableBuilder::finish`
    // keeps the first record of each endpoint pair, so file bonds pushed first
    // win over polymer links, which win over inferred grid bonds.
    let mut output = BondTableBuilder::new();
    for bond in structure.data().bonds.iter() {
        output.push(bond);
    }

    add_polymer_links(structure, &mut output);
    let thresholds: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| atom.element().and_then(element_threshold))
        .collect();
    let largest = thresholds.iter().flatten().copied().fold(0.0_f32, f32::max);
    let maximum = maximum_bonding_distance(largest, largest);
    if maximum > 0.0 && maximum.is_finite() {
        add_grid_bonds(structure, &thresholds, maximum, &mut output, context)?;
    }

    let mut data = structure.data().clone();
    data.bonds = output.finish();
    Ok(Structure::new(data))
}

/// Whether two atoms at a squared distance in ångström² form a plausible
/// covalent pair under the default bond-perception distance window.
///
/// Alternate conformers and H–H pairs are excluded. Metal coordination is not
/// covalent connectivity and is excluded even though default perception can
/// retain coordination edges. This is a geometric perception, not bond-order
/// assignment or a substitute for deposited connectivity.
#[must_use]
pub fn covalent_pair(left: AtomRef<'_>, right: AtomRef<'_>, distance_squared: f32) -> bool {
    let (Some(left_element), Some(right_element)) = (left.element(), right.element()) else {
        return false;
    };
    if is_metal(left_element.symbol()) || is_metal(right_element.symbol()) {
        return false;
    }
    let (Some(first), Some(second)) = (
        element_threshold(left_element),
        element_threshold(right_element),
    ) else {
        return false;
    };
    in_bonding_window(first, second, distance_squared) && pair_compatible(left, right)
}

/// One atom that can bond, stored where the grid keeps it.
#[derive(Clone, Copy)]
struct GridAtom {
    atom: u32,
    position: [f32; 3],
    threshold: f32,
}

/// Fewer grid atoms than this are searched on one block: a thread is not worth
/// waking for a structure that small.
const PARALLEL_GRID_ATOMS: usize = 32_768;

/// Grid cells searched per block of the parallel pass.
const CELLS_PER_BLOCK: usize = 4096;

fn add_grid_bonds(
    structure: &Structure,
    thresholds: &[Option<f32>],
    cutoff: f32,
    output: &mut BondTableBuilder,
    context: &ExecutionContext,
) -> Result<(), Diagnostic> {
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
    let single_block = entries.len() < PARALLEL_GRID_ATOMS;
    let grid = CellGrid::build(entries);
    let cells = grid.cell_count();
    let plan = if single_block {
        BlockPlan::new(cells, cells.max(1))
    } else {
        BlockPlan::new(cells, CELLS_PER_BLOCK)
    };
    // Blocks cover fixed cell ranges and are appended in block order, so the
    // records reach the builder in the order a serial walk would produce.
    let blocks = map_blocks_in(plan, context, |_, range| {
        let mut records = Vec::new();
        grid.visit_cells(range, |own, neighbourhood| {
            search_cell(structure, &grid, own, neighbourhood, &mut records);
        });
        records
    })
    .map_err(|_| {
        Diagnostic::new(Code::E1901).with_message("a worker thread panicked during bond perception")
    })?;
    for record in blocks.into_iter().flatten() {
        output.push(record);
    }
    Ok(())
}

/// Appends the bonds between one cell's atoms and everything in its neighbourhood.
fn search_cell(
    structure: &Structure,
    grid: &CellGrid<GridAtom>,
    own: Range<usize>,
    neighbourhood: &[Range<usize>],
    records: &mut Vec<BondRecord>,
) {
    for first in &grid.items()[own] {
        for range in neighbourhood {
            for second in &grid.items()[range.clone()] {
                if second.atom <= first.atom {
                    continue;
                }
                // Resolve thresholds once per atom and reject distant pairs
                // before materialising handles. covalent_pair shares this window.
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
                // Perception retains coordination; only covalent_pair excludes
                // metals. Both share conformer and H–H eligibility.
                if !pair_compatible(left, right) {
                    continue;
                }
                records.push(BondRecord {
                    atom_a: AtomIndex::new(first.atom),
                    atom_b: AtomIndex::new(second.atom),
                    order: intra_residue_order(left, right),
                    provenance: BondProvenance::InferredDistance,
                });
            }
        }
    }
}

fn add_polymer_links(structure: &Structure, output: &mut BondTableBuilder) {
    for chain in structure.data().chains() {
        let residues: Vec<_> = chain.residues().collect();
        for pair in residues.windows(2) {
            let [left_residue, right_residue] = pair else {
                continue;
            };
            if let Some((left, right)) = left_residue.atom("C").zip(right_residue.atom("N"))
                && within(left, right, 1.75)
            {
                add_link(left, right, output);
            }
            let left = ["O3'", "O3*", "O3"]
                .iter()
                .find_map(|name| left_residue.atom(name));
            if let Some((left, right)) = left.zip(right_residue.atom("P"))
                && within(left, right, 1.9)
            {
                add_link(left, right, output);
            }
        }
    }
}

fn add_link(left: AtomRef<'_>, right: AtomRef<'_>, output: &mut BondTableBuilder) {
    output.push(BondRecord {
        atom_a: left.index(),
        atom_b: right.index(),
        order: BondOrder::Polymeric,
        provenance: BondProvenance::InferredDistance,
    });
}

/// Whether a distance lies between the closest plausible contact and the sum of
/// the two elements' bonding thresholds.
fn in_bonding_window(left_threshold: f32, right_threshold: f32, distance_squared: f32) -> bool {
    let maximum = maximum_bonding_distance(left_threshold, right_threshold);
    distance_squared >= 0.16 && distance_squared <= maximum * maximum
}

fn maximum_bonding_distance(left_threshold: f32, right_threshold: f32) -> f32 {
    (left_threshold + right_threshold) / 1.95
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
