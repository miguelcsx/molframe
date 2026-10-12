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
use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::optional::OptionalSymbol;
use molframe_core::parallel::{
    BlockExecutionError, BlockPlan, map_blocks_in, try_for_each_block_in,
};
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use molframe_core::{AltId, SymbolId};
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
    let prepared = prepare_atoms(structure, context)?;
    let largest = prepared
        .iter()
        .filter_map(|atom| atom.threshold)
        .fold(0.0_f32, f32::max);
    let maximum = maximum_bonding_distance(largest, largest);
    if maximum > 0.0 && maximum.is_finite() {
        add_grid_bonds(structure, &prepared, maximum, &mut output, context)?;
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

/// Bonds an atom is charged for when sizing a block; a plan, not a limit.
const EXPECTED_BONDS_PER_ATOM: usize = 8;

/// Grid cells searched per block of the parallel pass.
const CELLS_PER_BLOCK: usize = 4096;

/// Atoms read per block while the per-atom facts are gathered.
const PREPARE_ATOMS_PER_BLOCK: usize = 16_384;

/// What the pair test needs from one atom, read once instead of once per pair.
///
/// Reading these through an atom handle costs a chunk lookup per field, which
/// dominated the search when it was done for every candidate pair.
#[derive(Clone, Copy, Default)]
struct AtomFacts {
    threshold: Option<f32>,
    hydrogen: bool,
    alt: AltId,
    residue: Option<ResidueIndex>,
    name: Option<SymbolId>,
    component: Option<SymbolId>,
}

/// Reads [`AtomFacts`] for every atom, in atom order, on the worker budget.
fn prepare_atoms(
    structure: &Structure,
    context: &ExecutionContext,
) -> Result<Vec<AtomFacts>, Diagnostic> {
    let count = structure.atom_count() as usize;
    let plan = BlockPlan::new(count, PREPARE_ATOMS_PER_BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range: Range<usize>| {
        range
            .map(|index| {
                u32::try_from(index)
                    .ok()
                    .and_then(|index| structure.atom(AtomIndex::new(index)))
                    .map_or_else(AtomFacts::default, |atom| AtomFacts {
                        threshold: atom.element().and_then(element_threshold),
                        hydrogen: atom
                            .element()
                            .is_some_and(molframe_core::Element::is_hydrogen),
                        alt: atom.alt_id().unwrap_or(AltId::BLANK),
                        residue: atom.residue().map(ResidueRef::index),
                        name: atom.name_symbol(),
                        component: atom.component_id(),
                    })
            })
            .collect::<Vec<_>>()
    })
    .map_err(|_| {
        Diagnostic::new(Code::E1901)
            .with_message("a worker thread panicked during bond perception")
    })?;
    let mut facts = Vec::with_capacity(count);
    for block in blocks {
        facts.extend(block);
    }
    Ok(facts)
}

fn add_grid_bonds(
    structure: &Structure,
    facts: &[AtomFacts],
    cutoff: f32,
    output: &mut BondTableBuilder,
    context: &ExecutionContext,
) -> Result<(), Diagnostic> {
    let mut entries = Vec::with_capacity(facts.len());
    for (atom, (position, fact)) in structure.positions().iter().zip(facts).enumerate() {
        let (Some(threshold), Some(cell), Ok(atom)) =
            (fact.threshold, cell_for(*position, cutoff), u32::try_from(atom))
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
    // records reach the builder in the order a serial walk would produce. Each
    // block is charged for the bonds its atoms are expected to hold.
    let block_cells = if single_block {
        cells.max(1)
    } else {
        CELLS_PER_BLOCK
    };
    let atoms_per_block = expected_block_atoms(grid.items().len(), cells, block_cells);
    let bytes_per_block = atoms_per_block
        .saturating_mul(EXPECTED_BONDS_PER_ATOM)
        .saturating_mul(size_of::<BondRecord>());
    let searched = try_for_each_block_in(
        plan,
        context,
        bytes_per_block,
        |_, range| {
            let mut records = Vec::new();
            grid.visit_cells(range, |own, neighbourhood| {
                search_cell(structure, facts, &grid, own, neighbourhood, &mut records);
            });
            Ok::<_, Diagnostic>(records)
        },
        |records| {
            for record in records {
                output.push(record);
            }
            Ok(())
        },
    );
    match searched {
        Ok(()) => Ok(()),
        // Parallel blocks are a speed-up, not a requirement: when the budget
        // cannot hold their buffers the same cells are walked on this thread,
        // in the same order, holding only one cell's bonds at a time. Bonds
        // an earlier window already pushed repeat here; the builder keeps
        // the first record of each pair. Reading
        // a structure never fails for want of search scratch.
        Err(BlockExecutionError::Memory(_)) => {
            let mut records = Vec::new();
            grid.visit_cells(0..cells, |own, neighbourhood| {
                search_cell(structure, facts, &grid, own, neighbourhood, &mut records);
                for record in records.drain(..) {
                    output.push(record);
                }
            });
            if context.cancellation().is_cancelled() {
                return Err(Diagnostic::new(Code::E1904)
                    .with_message("bond perception was cancelled"));
            }
            Ok(())
        }
        Err(BlockExecutionError::Cancelled) => {
            Err(Diagnostic::new(Code::E1904).with_message("bond perception was cancelled"))
        }
        Err(BlockExecutionError::Worker(_)) => Err(Diagnostic::new(Code::E1901)
            .with_message("a worker thread panicked during bond perception")),
        Err(BlockExecutionError::Operation(error)) => Err(error),
    }
}

/// Atoms expected in one block: the block's share of the cells.
///
/// A dense grid counts the empty slots of its bounding box. Rounding atoms per
/// cell up before scaling charged one atom for every empty cell and refused a
/// sparse 12k-atom structure for 115 MB it never needed, so the share is
/// scaled first and rounded once.
fn expected_block_atoms(atoms: usize, cells: usize, block_cells: usize) -> usize {
    let cells = cells.max(1);
    atoms
        .saturating_mul(block_cells.min(cells))
        .div_ceil(cells)
}

/// Appends the bonds between one cell's atoms and everything in its neighbourhood.
fn search_cell(
    structure: &Structure,
    facts: &[AtomFacts],
    grid: &CellGrid<GridAtom>,
    own: Range<usize>,
    neighbourhood: &[Range<usize>],
    records: &mut Vec<BondRecord>,
) {
    for first in &grid.items()[own] {
        let Some(left) = facts.get(first.atom as usize) else {
            continue;
        };
        for range in neighbourhood {
            for second in &grid.items()[range.clone()] {
                if second.atom <= first.atom {
                    continue;
                }
                // Resolve thresholds once per atom and reject distant pairs
                // before reading anything else. covalent_pair shares this window.
                let distance_squared = squared_distance(first.position, second.position);
                if !in_bonding_window(first.threshold, second.threshold, distance_squared) {
                    continue;
                }
                let Some(right) = facts.get(second.atom as usize) else {
                    continue;
                };
                // Perception retains coordination; only covalent_pair excludes
                // metals. Both share conformer and H–H eligibility.
                if !facts_compatible(left, right) {
                    continue;
                }
                records.push(BondRecord {
                    atom_a: AtomIndex::new(first.atom),
                    atom_b: AtomIndex::new(second.atom),
                    order: facts_order(structure, left, right),
                    provenance: BondProvenance::InferredDistance,
                    // Distance inference says two atoms are bonded; it says nothing
                    // about a deposited connection class.
                    kind: OptionalSymbol::NONE,
                });
            }
        }
    }
}

/// [`pair_compatible`] on facts read ahead of the search.
fn facts_compatible(left: &AtomFacts, right: &AtomFacts) -> bool {
    !(left.hydrogen && right.hydrogen)
        && (left.alt.is_blank() || right.alt.is_blank() || left.alt == right.alt)
}

/// [`intra_residue_order`] on facts read ahead of the search.
fn facts_order(structure: &Structure, left: &AtomFacts, right: &AtomFacts) -> BondOrder {
    let (Some(left_residue), Some(right_residue)) = (left.residue, right.residue) else {
        return BondOrder::Single;
    };
    if left_residue != right_residue {
        return BondOrder::Single;
    }
    let dictionary = &structure.data().dictionary;
    let (Some(component), Some(left_name), Some(right_name)) = (
        left.component.and_then(|symbol| dictionary.resolve(symbol)),
        left.name.and_then(|symbol| dictionary.resolve(symbol)),
        right.name.and_then(|symbol| dictionary.resolve(symbol)),
    ) else {
        return BondOrder::Single;
    };
    standard_bond_order(component, left_name, right_name)
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
        // The linkage is a polymer connection the distance implies, not a
        // deposited `struct_conn` row with a class of its own.
        kind: OptionalSymbol::NONE,
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

#[cfg(test)]
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

/// The bonding reach of one element, in ångström, or `None` when this reader
/// has no reference value for it.
///
/// An element with no entry takes part in no inferred bond at all: perception
/// leaves it to deposited connectivity rather than guessing a radius. Selenium
/// is the case that matters in practice — selenomethionine's `SE` is a real
/// covalent substituent — but the same list governs every other element, and
/// widening it is a data decision, not a code one.
fn element_threshold(element: molframe_core::Element) -> Option<f32> {
    match element.symbol() {
        "H" => Some(1.42),
        "C" => Some(1.75),
        "N" => Some(1.60),
        "O" => Some(1.52),
        "P" => Some(2.00),
        "S" => Some(1.90),
        "F" | "Cl" | "Br" | "I" => Some(1.80),
        symbol if is_metal(symbol) => metal_threshold(symbol),
        _ => None,
    }
}

/// The coordination reach of one metal, in ångström, on the same scale as the
/// non-metal thresholds.
///
/// A metal–ligand bond is accepted up to `(metal + ligand) / 1.95`, so a metal's
/// value is chosen so that the longest ordinary coordination distance of that
/// metal to a nitrogen, oxygen or sulfur donor still falls inside the window. One
/// shared value cannot do it: a heme iron sits 1.9–2.3 Å from its pyrrole
/// nitrogens and a calcium 2.3–2.6 Å from its oxygens, so a single compromise
/// either drops iron–nitrogen bonds that are really there or lets small metals
/// reach through to neighbours that are not coordinated.
fn metal_threshold(symbol: &str) -> Option<f32> {
    match symbol {
        "Li" => Some(2.80),
        "Na" => Some(3.30),
        "K" => Some(3.90),
        "Mg" | "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Au" | "Cr" | "V" => Some(2.90),
        "Ca" => Some(3.30),
        "Zn" => Some(2.85),
        "Cd" | "Mo" | "W" | "Ag" => Some(3.10),
        "Hg" => Some(3.00),
        "Al" => Some(2.80),
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
            | "Cr"
            | "V"
            | "W"
    )
}

#[cfg(test)]
#[path = "bonds_tests.rs"]
mod tests;
