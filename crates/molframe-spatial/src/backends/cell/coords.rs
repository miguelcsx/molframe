//! Cell-coordinate arithmetic and member lookup shared by every grid walk.
//!
//! Every function is `O(1)` except the sparse lookup, which binary-searches the
//! occupied-cell list.

use super::{CellId, Grid};
use crate::SpatialError;

/// Converts a Cartesian point into clamped cell coordinates.
///
/// Runtime and auxiliary space are `O(1)`.
pub(super) fn raw_cell(
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    point: [f32; 3],
) -> Result<[usize; 3], SpatialError> {
    let convert = |axis: usize| {
        // Differences are formed in double precision: at a large extent the
        // single-precision subtraction would round cluster atoms into wrong cells.
        let scaled = (f64::from(point[axis]) - f64::from(origin[axis])) / f64::from(edge);
        if scaled.is_finite() && scaled <= 0.0 {
            return Ok(0);
        }
        let last = dims[axis]
            .checked_sub(1)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        crate::numeric::floor_usize(scaled)
            .map(|raw| raw.min(last))
            .ok_or(SpatialError::NumericRangeExceeded)
    };
    Ok([convert(0)?, convert(1)?, convert(2)?])
}

/// Converts a Cartesian point into this grid's cell coordinate.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(super) fn cell_of(grid: &Grid, point: [f32; 3]) -> Result<[usize; 3], SpatialError> {
    raw_cell(grid.origin, grid.edge, grid.dims, point)
}

/// Converts three-dimensional cell coordinates to a linear cell index.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(super) fn linear(dims: [usize; 3], cell: [usize; 3]) -> CellId {
    let wide = |value: usize| value as CellId;
    (wide(cell[0]) * wide(dims[1]) + wide(cell[1])) * wide(dims[2]) + wide(cell[2])
}

/// [`linear`] as a position in a dense table, which `choose_layout` bounded.
#[inline]
pub(super) fn linear_index(dims: [usize; 3], cell: [usize; 3]) -> Result<usize, SpatialError> {
    usize::try_from(linear(dims, cell)).map_err(|_| SpatialError::NumericRangeExceeded)
}

/// Returns a neighbouring linear cell index when it remains in bounds.
///
/// Runtime and auxiliary space are `O(1)`.
pub(super) fn neighbour_cell(
    dims: [usize; 3],
    centre: [usize; 3],
    delta: [isize; 3],
) -> Option<CellId> {
    let mut cell = [0usize; 3];

    for axis in 0..3 {
        let coordinate = centre[axis].checked_add_signed(delta[axis])?;

        if coordinate >= dims[axis] {
            return None;
        }

        cell[axis] = coordinate;
    }

    Some(linear(dims, cell))
}

/// Returns the contiguous target-member slice belonging to `cell`.
///
/// Invalid internal offsets defensively produce an empty slice.
pub(super) fn members_of_cell(grid: &Grid, cell: CellId) -> Option<&[u32]> {
    if grid.sparse {
        let Ok(slot) = grid.non_empty_cells.binary_search(&cell) else {
            return Some(&[]);
        };
        let &start = grid.offsets.get(slot)?;
        let &end = grid.offsets.get(slot + 1)?;
        return grid.members.get(start..end);
    }
    let cell = usize::try_from(cell).ok()?;
    let &start = grid.offsets.get(cell)?;
    let next = cell.checked_add(1)?;
    let &end = grid.offsets.get(next)?;

    grid.members.get(start..end)
}
