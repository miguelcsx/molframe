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

/// One atom that can bond, stored where the grid keeps it.
#[derive(Clone, Copy)]
struct GridAtom {
    atom: u32,
    position: [f32; 3],
    threshold: f32,
}

/// The most cells a dense grid may hold; a sparser, larger span is searched
/// through its occupied cells instead.
const DENSE_CELL_LIMIT: usize = 1 << 22;

/// Bondable atoms ordered by grid cell.
///
/// Positions and thresholds are copied into cell order, so scanning a
/// neighbouring cell reads one contiguous run instead of jumping through the
/// structure. An atom without a threshold or a finite position can never bond
/// and is left out.
struct CellGrid {
    atoms: Vec<GridAtom>,
    layout: Layout,
}

/// How the cells of a [`CellGrid`] are found.
enum Layout {
    /// One slot per cell of the bounding box: a neighbour is an index away.
    Dense { offsets: Vec<u32>, dims: [usize; 3] },
    /// Only occupied cells, searched by key, for a span too large to tabulate.
    Sparse {
        cells: Vec<(CellKey, std::ops::Range<usize>)>,
    },
}

impl CellGrid {
    fn build(positions: &[[f32; 3]], thresholds: &[Option<f32>], cutoff: f32) -> Self {
        let mut entries = Vec::with_capacity(positions.len());
        for (atom, (position, threshold)) in positions.iter().zip(thresholds).enumerate() {
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
        match Self::dense(&entries) {
            Some(grid) => grid,
            None => Self::sparse(entries),
        }
    }

    /// Counting-sorts the entries into the cells of their bounding box.
    ///
    /// Returns `None` when the box would hold more cells than the limit.
    fn dense(entries: &[(CellKey, GridAtom)]) -> Option<Self> {
        let first = entries.first()?.0.0;
        let (mut low, mut high) = (first, first);
        for (cell, _) in entries {
            for axis in 0..3 {
                low[axis] = low[axis].min(cell.0[axis]);
                high[axis] = high[axis].max(cell.0[axis]);
            }
        }
        let mut dims = [0_usize; 3];
        let mut total = 1_usize;
        for axis in 0..3 {
            let span = i64::from(high[axis]) - i64::from(low[axis]) + 1;
            dims[axis] = usize::try_from(span).ok()?;
            total = total.checked_mul(dims[axis])?;
        }
        if total > DENSE_CELL_LIMIT {
            return None;
        }
        let slot = |cell: &CellKey| -> usize {
            // Every entry lies inside the bounding box, so the offsets from its
            // corner are never negative; a failed conversion cannot occur.
            let at =
                |axis: usize| match usize::try_from(i64::from(cell.0[axis]) - i64::from(low[axis]))
                {
                    Ok(offset) => offset,
                    Err(_) => 0,
                };
            (at(0) * dims[1] + at(1)) * dims[2] + at(2)
        };
        let mut offsets = vec![0_u32; total + 1];
        for (cell, _) in entries {
            offsets[slot(cell) + 1] += 1;
        }
        for index in 1..offsets.len() {
            offsets[index] += offsets[index - 1];
        }
        // Filling in entry order keeps each cell's atoms ascending.
        let mut fill: Vec<u32> = offsets[..total].to_vec();
        let mut atoms = vec![entries.first()?.1; entries.len()];
        for (cell, atom) in entries {
            let target = &mut fill[slot(cell)];
            atoms[*target as usize] = *atom;
            *target += 1;
        }
        Some(Self {
            atoms,
            layout: Layout::Dense { offsets, dims },
        })
    }

    fn sparse(mut entries: Vec<(CellKey, GridAtom)>) -> Self {
        entries.sort_unstable_by_key(|(cell, atom)| (*cell, atom.atom));
        let mut cells: Vec<(CellKey, std::ops::Range<usize>)> = Vec::new();
        for (position, (cell, _)) in entries.iter().enumerate() {
            match cells.last_mut() {
                Some((last, range)) if last == cell => range.end = position + 1,
                _ => cells.push((*cell, position..position + 1)),
            }
        }
        Self {
            atoms: entries.into_iter().map(|(_, atom)| atom).collect(),
            layout: Layout::Sparse { cells },
        }
    }

    /// Calls `visit` with each occupied cell's atoms and the atoms of every
    /// occupied cell around it, itself included.
    fn for_each_cell(
        &self,
        mut visit: impl FnMut(std::ops::Range<usize>, &[std::ops::Range<usize>]),
    ) {
        let mut neighbourhood = Vec::with_capacity(27);
        match &self.layout {
            Layout::Dense { offsets, dims } => {
                let range_at = |x: usize, y: usize, z: usize| {
                    let slot = (x * dims[1] + y) * dims[2] + z;
                    offsets[slot] as usize..offsets[slot + 1] as usize
                };
                for x in 0..dims[0] {
                    for y in 0..dims[1] {
                        for z in 0..dims[2] {
                            let own = range_at(x, y, z);
                            if own.is_empty() {
                                continue;
                            }
                            neighbourhood.clear();
                            for nx in x.saturating_sub(1)..=(x + 1).min(dims[0] - 1) {
                                for ny in y.saturating_sub(1)..=(y + 1).min(dims[1] - 1) {
                                    for nz in z.saturating_sub(1)..=(z + 1).min(dims[2] - 1) {
                                        let range = range_at(nx, ny, nz);
                                        if !range.is_empty() {
                                            neighbourhood.push(range);
                                        }
                                    }
                                }
                            }
                            visit(own, &neighbourhood);
                        }
                    }
                }
            }
            Layout::Sparse { cells } => {
                for (base, own) in cells {
                    neighbourhood.clear();
                    for dx in -1..=1 {
                        for dy in -1..=1 {
                            for dz in -1..=1 {
                                let key = CellKey([
                                    base.0[0].saturating_add(dx),
                                    base.0[1].saturating_add(dy),
                                    base.0[2].saturating_add(dz),
                                ]);
                                if let Ok(index) = cells.binary_search_by_key(&key, |(key, _)| *key)
                                {
                                    neighbourhood.push(cells[index].1.clone());
                                }
                            }
                        }
                    }
                    visit(own.clone(), &neighbourhood);
                }
            }
        }
    }
}

fn add_grid_bonds(
    structure: &Structure,
    thresholds: &[Option<f32>],
    cutoff: f32,
    output: &mut BondTableBuilder,
    existing: &mut HashSet<(u32, u32)>,
) {
    let grid = CellGrid::build(structure.positions(), thresholds, cutoff);
    grid.for_each_cell(|own, neighbourhood| {
        for first in &grid.atoms[own] {
            for range in neighbourhood {
                for second in &grid.atoms[range.clone()] {
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
