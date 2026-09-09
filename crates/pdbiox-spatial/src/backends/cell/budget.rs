//! Reserve construction peak before allocating a cell index.

use super::{CellGrid, CellList, finite_bounds, grid_geometry, periodic_grid};
use crate::{CellGridOptions, PeriodicBox, SpatialError};
use pdbiox_core::ExecutionContext;

impl<'a> CellList<'a> {
    /// Builds an index whose construction peak and retained buffers share the
    /// execution account. Dropping the index releases its retained charge.
    ///
    /// # Errors
    ///
    /// Returns invalid geometry, cancellation or memory admission errors before
    /// allocating the grid buffers.
    pub fn build_in(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        cutoff: f32,
        periodic: Option<&'a PeriodicBox>,
        options: CellGridOptions,
        context: &ExecutionContext,
    ) -> Result<Self, SpatialError> {
        if context.cancellation().is_cancelled() {
            return Err(SpatialError::Cancelled);
        }
        super::validate_cutoff(cutoff)?;
        super::validate_indices(targets, positions.len())?;
        options.validate()?;
        let bytes = workspace_bytes(positions, targets, cutoff, periodic, options)?;
        let mut reservation = context.try_reserve(bytes)?;
        let mut index = Self::build_with_options(positions, targets, cutoff, periodic, options)?;
        let retained = match &index.grid {
            Some(CellGrid::Cartesian(grid)) => grid.retained_bytes(),
            Some(CellGrid::Periodic(grid)) => grid.retained_bytes(),
            None => 0,
        };
        reservation.shrink_to(retained);
        index.reservation = Some(reservation);
        Ok(index)
    }
}

fn workspace_bytes(
    positions: &[[f32; 3]],
    targets: &[u32],
    cutoff: f32,
    periodic: Option<&PeriodicBox>,
    options: CellGridOptions,
) -> Result<usize, SpatialError> {
    let (cells, atoms, occupied) = if let Some(periodic) = periodic {
        let dims = periodic_grid::periodic_dimensions(periodic, cutoff, options)?;
        (periodic_grid::cell_count(dims)?, targets.len(), 0)
    } else if let Some(bounds) = finite_bounds(positions, targets) {
        let cells = if cutoff == 0.0 {
            1
        } else {
            grid_geometry(bounds.min, bounds.max, cutoff, options)?.2
        };
        (cells, bounds.finite_count, cells.min(bounds.finite_count))
    } else {
        return Ok(0);
    };
    // Counts and prefix offsets coexist with members and occupied-cell indices.
    // Each vector reserves its final capacity, so construction never doubles it.
    cells
        .checked_mul(2)
        .and_then(|count| count.checked_add(1))
        .and_then(|count| count.checked_add(occupied))
        .and_then(|count| count.checked_mul(size_of::<usize>()))
        .and_then(|bytes| {
            atoms
                .checked_mul(size_of::<u32>())
                .and_then(|members| bytes.checked_add(members))
        })
        .ok_or(SpatialError::NumericRangeExceeded)
}

#[cfg(test)]
#[path = "budget_tests.rs"]
mod tests;

impl super::Grid {
    fn retained_bytes(&self) -> usize {
        (self.offsets.capacity() + self.non_empty_cells.capacity()) * size_of::<usize>()
            + self.members.capacity() * size_of::<u32>()
    }
}
