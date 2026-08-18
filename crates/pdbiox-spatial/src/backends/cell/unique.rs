//! Pair enumeration specialised for identical query and target selections.

use super::{CellGrid, CellList, Grid, NEIGHBOUR_OFFSETS};
use crate::brute::{canonicalise, distance_squared};
use crate::{NeighborPair, SpatialError};

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
        let members =
            super::members_of_cell(grid, cell).ok_or(SpatialError::NumericRangeExceeded)?;
        append_same_cell_pairs(positions, members, cutoff_squared, emit)?;

        let centre = cell_coordinates(grid.dims, cell).ok_or(SpatialError::NumericRangeExceeded)?;
        for delta in NEIGHBOUR_OFFSETS {
            let Some(neighbour) = super::neighbour_cell(grid.dims, centre, delta) else {
                continue;
            };
            if neighbour <= cell {
                continue;
            }
            let targets = super::members_of_cell(grid, neighbour)
                .ok_or(SpatialError::NumericRangeExceeded)?;
            append_cross_cell_pairs(positions, members, targets, cutoff_squared, emit)?;
        }
    }
    Ok(())
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
