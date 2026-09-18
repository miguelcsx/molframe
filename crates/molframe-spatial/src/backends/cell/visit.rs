//! Streaming cross-selection queries over a Cartesian cell grid.

use super::{Grid, NEIGHBOUR_OFFSETS, cell_of, members_of_cell, neighbour_cell};
use crate::SpatialError;
use crate::brute::finite;

/// Evaluates a cell-grid query and emits matches as they are found.
///
/// Each finite query atom visits at most 27 cells. Pairwise work inside those
/// cells depends on local spatial density.
pub(super) fn grid_for_each_pair<const UPPER: bool>(
    positions: &[[f32; 3]],
    grid: &Grid,
    query: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) -> Result<(), SpatialError> {
    for &atom in query {
        let Ok(index) = usize::try_from(atom) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        if !finite(position) {
            continue;
        }

        append_grid_pairs::<UPPER>(positions, grid, atom, position, cutoff_squared, emit)?;
    }
    Ok(())
}

/// Searches the 27 neighbouring cells for one query atom.
fn append_grid_pairs<const UPPER: bool>(
    positions: &[[f32; 3]],
    grid: &Grid,
    atom: u32,
    position: [f32; 3],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) -> Result<(), SpatialError> {
    let centre = cell_of(grid, position)?;

    for delta in NEIGHBOUR_OFFSETS {
        let Some(cell) = neighbour_cell(grid.dims, centre, delta) else {
            continue;
        };

        let members = members_of_cell(grid, cell).ok_or(SpatialError::NumericRangeExceeded)?;
        super::super::brute_simd::append_pairs::<UPPER>(
            positions,
            atom,
            position,
            members,
            cutoff_squared,
            emit,
        );
    }
    Ok(())
}
