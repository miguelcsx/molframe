//! Canonical Cartesian cell ownership for bounded block consumers.

use super::{CellGrid, CellList};
use crate::{NeighborPair, SpatialError};
use std::ops::Range;

impl CellList<'_> {
    pub(crate) fn owned_cell_count(&self) -> Option<usize> {
        match &self.grid {
            Some(CellGrid::Cartesian(grid)) if grid.largest_cell <= 256 => {
                Some(grid.non_empty_cells.len())
            }
            _ => None,
        }
    }

    // Query tiles split dense cells while rejecting reverse orientations before
    // distance arithmetic. Tile boundaries do not depend on worker count.
    pub(crate) fn for_each_upper_candidate(
        &self,
        query: &[u32],
        cutoff: f32,
        mut emit: impl FnMut(NeighborPair),
    ) -> Result<(), SpatialError> {
        super::validate_query_cutoff(cutoff, self.cutoff)?;
        super::validate_indices(query, self.positions.len())?;
        match &self.grid {
            Some(CellGrid::Cartesian(grid)) => super::visit::grid_for_each_pair::<true>(
                self.positions,
                grid,
                query,
                cutoff * cutoff,
                &mut |left, right, squared| emit(NeighborPair::new(left, right, squared)),
            ),
            Some(CellGrid::Periodic(grid)) => {
                let periodic = self.periodic.ok_or(SpatialError::InvalidCell)?;
                grid.for_each_pairs_same_selection(
                    self.positions,
                    query,
                    cutoff * cutoff,
                    periodic,
                    &mut emit,
                )
            }
            None => Ok(()),
        }
    }

    pub(crate) fn for_each_owned_cell(
        &self,
        cells: Range<usize>,
        cutoff: f32,
        mut emit: impl FnMut(NeighborPair),
    ) -> Result<(), SpatialError> {
        super::validate_query_cutoff(cutoff, self.cutoff)?;
        let Some(CellGrid::Cartesian(grid)) = &self.grid else {
            return Err(SpatialError::InvalidCell);
        };
        let cells = grid
            .non_empty_cells
            .get(cells)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        for &cell in cells {
            super::unique::emit_cell_pairs(self.positions, grid, cell, cutoff * cutoff, &mut emit)?;
        }
        Ok(())
    }
}
