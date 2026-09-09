//! Pair enumeration specialised for identical query and target selections.

use super::{CellGrid, CellList, Grid, NEIGHBOUR_OFFSETS};
use crate::brute::{canonicalise, distance_squared};
use crate::{NeighborPair, SpatialError};
use pdbiox_core::ExecutionContext;
use pdbiox_core::parallel::{BlockPlan, map_blocks_in};

/// Finds each unordered pair from one selection once.
///
/// The ordinary cell query visits both `(a, b)` and `(b, a)` when the query
/// and indexed target selections are identical. Restricting emission to
/// `a < b` avoids that duplicate distance work while retaining the public
/// sorted-and-unique result contract after canonicalisation.
pub(crate) fn pairs_same_selection(
    list: &CellList<'_>,
    query: &[u32],
    cutoff: f32,
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_same_selection_with_order(list, query, cutoff, true)
}

/// Finds each unordered pair once without sorting the result.
pub(crate) fn pairs_same_selection_unordered(
    list: &CellList<'_>,
    query: &[u32],
    cutoff: f32,
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_same_selection_with_order(list, query, cutoff, false)
}

/// Visits each unordered same-selection pair once without materialising it.
pub(crate) fn for_each_pairs_same_selection_unordered<F>(
    list: &CellList<'_>,
    query: &[u32],
    cutoff: f32,
    mut emit: F,
) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    super::validate_query_cutoff(cutoff, list.cutoff)?;
    super::validate_indices(query, list.positions.len())?;

    match &list.grid {
        Some(CellGrid::Cartesian(grid)) => {
            grid_for_each_pairs(list.positions, grid, cutoff * cutoff, &mut emit)
        }
        Some(CellGrid::Periodic(grid)) => {
            let Some(periodic) = list.periodic else {
                return Err(SpatialError::InvalidCell);
            };
            grid.for_each_pairs_same_selection(
                list.positions,
                query,
                cutoff * cutoff,
                periodic,
                &mut emit,
            )
        }
        None => Ok(()),
    }
}

fn pairs_same_selection_with_order(
    list: &CellList<'_>,
    query: &[u32],
    cutoff: f32,
    sort_result: bool,
) -> Result<Vec<NeighborPair>, SpatialError> {
    super::validate_query_cutoff(cutoff, list.cutoff)?;
    super::validate_indices(query, list.positions.len())?;

    match &list.grid {
        Some(CellGrid::Cartesian(grid)) => {
            let mut found = Vec::with_capacity(query.len());
            grid_for_each_pairs(list.positions, grid, cutoff * cutoff, &mut |pair| {
                found.push(pair);
            })?;
            if sort_result {
                canonicalise(&mut found);
            }
            Ok(found)
        }
        Some(CellGrid::Periodic(grid)) => {
            let Some(periodic) = list.periodic else {
                return Err(SpatialError::InvalidCell);
            };
            grid.pairs_same_selection(
                list.positions,
                query,
                cutoff * cutoff,
                periodic,
                sort_result,
            )
        }
        None => Ok(Vec::new()),
    }
}

/// Evaluates a Cartesian cell query while emitting only `atom < target`.
fn grid_for_each_pairs<F>(
    positions: &[[f32; 3]],
    grid: &Grid,
    cutoff_squared: f32,
    emit: &mut F,
) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    for &cell in &grid.non_empty_cells {
        emit_cell_pairs(positions, grid, cell, cutoff_squared, emit)?;
    }
    Ok(())
}

/// Emits every unordered pair owned by one non-empty cell: pairs inside the
/// cell, and cross-cell pairs with each higher-indexed neighbour. Each pair is
/// owned by exactly one cell, so partitioning `non_empty_cells` across threads
/// yields a disjoint emission with no double counting.
pub(super) fn emit_cell_pairs<F>(
    positions: &[[f32; 3]],
    grid: &Grid,
    cell: usize,
    cutoff_squared: f32,
    emit: &mut F,
) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    let members = super::members_of_cell(grid, cell).ok_or(SpatialError::NumericRangeExceeded)?;
    append_same_cell_pairs(positions, members, cutoff_squared, emit)?;

    let centre = cell_coordinates(grid.dims, cell).ok_or(SpatialError::NumericRangeExceeded)?;
    for delta in NEIGHBOUR_OFFSETS {
        let Some(neighbour) = super::neighbour_cell(grid.dims, centre, delta) else {
            continue;
        };
        if neighbour <= cell {
            continue;
        }
        let targets =
            super::members_of_cell(grid, neighbour).ok_or(SpatialError::NumericRangeExceeded)?;
        append_cross_cell_pairs(positions, members, targets, cutoff_squared, emit)?;
    }
    Ok(())
}

/// Finds each unordered same-selection pair across `workers` scoped threads.
///
/// `non_empty_cells` is divided by a block plan whose boundaries come from the
/// cell count, not the worker count; each cell owns a disjoint set of pairs, so
/// the concatenated result is exactly the serial set. `canonicalise` then sorts
/// it, making the output identical to [`pairs_same_selection`] for any worker
/// count. Non-Cartesian grids fall back to the serial sorted path.
pub(crate) fn pairs_same_selection_parallel(
    list: &CellList<'_>,
    query: &[u32],
    cutoff: f32,
    context: &ExecutionContext,
    sort_result: bool,
) -> Result<Vec<NeighborPair>, SpatialError> {
    /// Cells per block. Occupancy varies widely between a dense protein core
    /// and a solvent shell, so blocks are much smaller than one worker's share.
    const BLOCK_CELLS: usize = 64;
    /// Below this size, dispatch and per-block vectors cost more than the
    /// available cell parallelism on the calibrated Apple M4 workload.
    const MIN_PARALLEL_ATOMS: usize = 2_048;

    super::validate_query_cutoff(cutoff, list.cutoff)?;
    super::validate_indices(query, list.positions.len())?;

    let Some(CellGrid::Cartesian(grid)) = &list.grid else {
        return pairs_same_selection_with_order(list, query, cutoff, sort_result);
    };
    if query.len() < MIN_PARALLEL_ATOMS {
        return pairs_same_selection_with_order(list, query, cutoff, sort_result);
    }

    let cutoff_squared = cutoff * cutoff;
    let cells = &grid.non_empty_cells;
    let positions = list.positions;

    let plan = BlockPlan::new(cells.len(), BLOCK_CELLS);
    let produced = map_blocks_in(plan, context, |_, range| match cells.get(range) {
        Some(chunk) => collect_cells(positions, grid, chunk, cutoff_squared),
        None => Ok(Vec::new()),
    })
    .map_err(|_| SpatialError::WorkerPanicked)?;

    let mut found = Vec::new();
    for part in produced {
        found.extend(part?);
    }

    if sort_result {
        canonicalise(&mut found);
    }
    Ok(found)
}

/// Folds every same-selection pair into per-block accumulators.
///
/// Identical block partitioning to [`pairs_same_selection_parallel`], but each
/// block reduces its pairs as it finds them instead of collecting them. Peak
/// retained bytes then track the accumulator, not the pair count — the
/// difference between a reduction that fits in cache and one that needs
/// hundreds of gigabytes at a billion atoms (FR-516).
///
/// Accumulators are returned in block order, so merging them in the order given
/// produces the same result at any worker count (FR-515).
pub(crate) fn reduce_pairs_same_selection_parallel<T, I, F>(
    list: &CellList<'_>,
    query: &[u32],
    cutoff: f32,
    context: &ExecutionContext,
    init: I,
    fold: F,
) -> Result<Option<Vec<T>>, SpatialError>
where
    T: Send,
    I: Fn() -> T + Sync,
    F: Fn(&mut T, NeighborPair) + Sync,
{
    /// Cells per block, matching the collecting kernel so both partition alike.
    const BLOCK_CELLS: usize = 64;

    super::validate_query_cutoff(cutoff, list.cutoff)?;
    super::validate_indices(query, list.positions.len())?;

    // Only a Cartesian grid owns disjoint pairs per cell; other grids have no
    // block decomposition, so the caller falls back to the serial visitor.
    let Some(CellGrid::Cartesian(grid)) = &list.grid else {
        return Ok(None);
    };

    let cutoff_squared = cutoff * cutoff;
    let cells = &grid.non_empty_cells;
    let positions = list.positions;

    let plan = BlockPlan::new(cells.len(), BLOCK_CELLS);
    let produced = map_blocks_in(plan, context, |_, range| {
        let mut accumulator = init();
        let Some(chunk) = cells.get(range) else {
            return Ok(accumulator);
        };
        for &cell in chunk {
            emit_cell_pairs(positions, grid, cell, cutoff_squared, &mut |pair| {
                fold(&mut accumulator, pair);
            })?;
        }
        Ok(accumulator)
    })
    .map_err(|_| SpatialError::WorkerPanicked)?;

    produced
        .into_iter()
        .collect::<Result<Vec<T>, _>>()
        .map(Some)
}

/// Collects every pair owned by a contiguous slice of non-empty cells.
fn collect_cells(
    positions: &[[f32; 3]],
    grid: &Grid,
    cells: &[usize],
    cutoff_squared: f32,
) -> Result<Vec<NeighborPair>, SpatialError> {
    let mut found = Vec::new();
    for &cell in cells {
        emit_cell_pairs(positions, grid, cell, cutoff_squared, &mut |pair| {
            found.push(pair);
        })?;
    }
    Ok(found)
}

/// Converts a linear cell index back to three-dimensional coordinates.
#[inline]
fn cell_coordinates(dims: [usize; 3], cell: usize) -> Option<[usize; 3]> {
    let plane = dims[1].checked_mul(dims[2])?;
    if plane == 0 || cell >= dims[0].checked_mul(plane)? {
        return None;
    }
    let first = cell / plane;
    let remainder = cell % plane;
    Some([first, remainder / dims[2], remainder % dims[2]])
}

/// Compares every pair inside one cell once by using canonical atom IDs.
fn append_same_cell_pairs<F>(
    positions: &[[f32; 3]],
    members: &[u32],
    cutoff_squared: f32,
    emit: &mut F,
) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    for &atom in members {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        for &target in members {
            if atom >= target {
                continue;
            }
            let target_index =
                usize::try_from(target).map_err(|_| SpatialError::NumericRangeExceeded)?;
            let Some(target_position) = positions.get(target_index).copied() else {
                continue;
            };

            let squared = distance_squared(position, target_position, None);

            if squared <= cutoff_squared {
                emit(NeighborPair::new(atom, target, squared));
            }
        }
    }
    Ok(())
}

/// Compares every atom pair between two distinct neighbouring cells once.
fn append_cross_cell_pairs<F>(
    positions: &[[f32; 3]],
    members: &[u32],
    targets: &[u32],
    cutoff_squared: f32,
    emit: &mut F,
) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    for &atom in members {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        for &target in targets {
            let target_index =
                usize::try_from(target).map_err(|_| SpatialError::NumericRangeExceeded)?;
            let Some(target_position) = positions.get(target_index).copied() else {
                continue;
            };

            let squared = distance_squared(position, target_position, None);

            if squared <= cutoff_squared {
                emit(NeighborPair::new(atom, target, squared));
            }
        }
    }
    Ok(())
}
