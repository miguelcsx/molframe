//! Contiguous fixed-radius cell list.

use crate::brute::{canonicalise, finite};
use crate::{CellGridOptions, NeighborPair, PeriodicBox, SpatialError};

#[path = "cell/budget.rs"]
mod budget;
#[path = "cell/build.rs"]
mod build;
use build::build_grid;
#[path = "cell/coords.rs"]
mod coords;
use coords::{cell_of, linear, linear_index, members_of_cell, neighbour_cell, raw_cell};
#[path = "cell/geometry.rs"]
mod geometry;
#[path = "cell/owned.rs"]
mod owned;
use geometry::{Layout, choose_layout, finite_bounds};
#[path = "cell/periodic.rs"]
mod periodic_grid;
use periodic_grid::PeriodicGrid;
#[path = "cell/unique.rs"]
mod unique;
pub(crate) use unique::{
    for_each_pairs_same_selection_unordered, pairs_same_selection, pairs_same_selection_parallel,
    pairs_same_selection_unordered, reduce_pairs_same_selection_parallel,
};
#[path = "cell/visit.rs"]
mod visit;
use visit::grid_for_each_pair;

const NEIGHBOUR_OFFSETS: [[isize; 3]; 27] = [
    [-1, -1, -1],
    [-1, -1, 0],
    [-1, -1, 1],
    [-1, 0, -1],
    [-1, 0, 0],
    [-1, 0, 1],
    [-1, 1, -1],
    [-1, 1, 0],
    [-1, 1, 1],
    [0, -1, -1],
    [0, -1, 0],
    [0, -1, 1],
    [0, 0, -1],
    [0, 0, 0],
    [0, 0, 1],
    [0, 1, -1],
    [0, 1, 0],
    [0, 1, 1],
    [1, -1, -1],
    [1, -1, 0],
    [1, -1, 1],
    [1, 0, -1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, -1],
    [1, 1, 0],
    [1, 1, 1],
];

/// Linear id of a cell. Wide enough for three axes of any single-precision
/// extent, so a sparse grid never has to coarsen its cells for a far outlier.
type CellId = u128;

#[derive(Debug)]
struct Grid {
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    offsets: Vec<usize>,
    members: Vec<u32>,
    /// Linear ids of occupied cells in ascending order.
    non_empty_cells: Vec<CellId>,
    largest_cell: usize,
    /// When set, `offsets` is indexed by position in `non_empty_cells` rather
    /// than by linear cell id, so memory is proportional to the atom count.
    sparse: bool,
}

#[derive(Debug)]
enum CellGrid {
    Cartesian(Grid),
    Periodic(PeriodicGrid),
}

/// Fixed-radius grid built in O(n) with one contiguous member array.
#[derive(Debug)]
pub struct CellList<'a> {
    positions: &'a [[f32; 3]],
    cutoff: f32,
    periodic: Option<PeriodicBox>,
    grid: Option<CellGrid>,
    reservation: Option<molframe_core::MemoryReservation>,
}

impl<'a> CellList<'a> {
    /// Builds an index over `targets`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid cutoff or target index.
    pub fn build(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        cutoff: f32,
        periodic: Option<PeriodicBox>,
    ) -> Result<Self, SpatialError> {
        Self::build_with_options(
            positions,
            targets,
            cutoff,
            periodic,
            CellGridOptions::default(),
        )
    }

    /// Builds an index under an explicit cell-grid memory policy.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid options, cutoff or target index.
    pub fn build_with_options(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        cutoff: f32,
        periodic: Option<PeriodicBox>,
        options: CellGridOptions,
    ) -> Result<Self, SpatialError> {
        validate_cutoff(cutoff)?;
        validate_indices(targets, positions.len())?;
        options.validate()?;

        let grid = match periodic.as_ref() {
            Some(periodic) => Some(CellGrid::Periodic(PeriodicGrid::build(
                positions, targets, cutoff, periodic, options,
            )?)),
            None => build_grid(positions, targets, cutoff, options)?.map(CellGrid::Cartesian),
        };

        Ok(Self {
            positions,
            cutoff,
            periodic,
            grid,
            reservation: None,
        })
    }

    /// Finds pairs for the indexed targets.
    ///
    /// # Errors
    ///
    /// Returns an error when the requested cutoff exceeds the build cutoff or
    /// a query index is out of range.
    pub fn pairs(&self, query: &[u32], cutoff: f32) -> Result<Vec<NeighborPair>, SpatialError> {
        let mut found = Vec::with_capacity(query.len());
        self.for_each_pair(query, cutoff, |pair| found.push(pair))?;
        canonicalise(&mut found);
        Ok(found)
    }

    /// Visits pairs for the indexed targets without retaining a result vector.
    pub(crate) fn for_each_pair(
        &self,
        query: &[u32],
        cutoff: f32,
        mut emit: impl FnMut(NeighborPair),
    ) -> Result<(), SpatialError> {
        self.for_each_candidate(query, cutoff, |atom, target, squared| {
            emit(NeighborPair::new(atom, target, squared));
        })
    }

    /// Visits oriented query-target candidates without retaining pairs.
    ///
    /// Self-pairs are omitted. Overlapping query and target sets may emit both
    /// orientations; callers computing symmetric results must choose ownership.
    ///
    /// # Errors
    ///
    /// Returns invalid query indices or a cutoff larger than the indexed cutoff.
    pub fn for_each_candidate(
        &self,
        query: &[u32],
        cutoff: f32,
        mut emit: impl FnMut(u32, u32, f32),
    ) -> Result<(), SpatialError> {
        validate_query_cutoff(cutoff, self.cutoff)?;
        validate_indices(query, self.positions.len())?;

        match &self.grid {
            Some(CellGrid::Cartesian(grid)) => {
                grid_for_each_pair::<false>(self.positions, grid, query, cutoff * cutoff, &mut emit)
            }
            Some(CellGrid::Periodic(grid)) => {
                let Some(periodic) = &self.periodic else {
                    return Err(SpatialError::InvalidCell);
                };
                grid.for_each_pair(self.positions, query, cutoff * cutoff, periodic, &mut emit)
            }
            None => Ok(()),
        }
    }

    /// Visits indexed atoms within `cutoff` of an arbitrary Cartesian point.
    ///
    /// Emits target indices and squared distances in ångström² without
    /// allocation. Unlike atom queries, a coincident target is included.
    /// Periodic indices use the same minimum-image metric as pair queries.
    ///
    /// # Errors
    ///
    /// Returns an error for a nonfinite point or a cutoff exceeding the index.
    pub fn for_each_neighbor(
        &self,
        position: [f32; 3],
        cutoff: f32,
        mut emit: impl FnMut(u32, f32),
    ) -> Result<(), SpatialError> {
        validate_query_cutoff(cutoff, self.cutoff)?;
        if !finite(position) {
            return Err(SpatialError::NumericRangeExceeded);
        }
        match &self.grid {
            Some(CellGrid::Cartesian(grid)) => visit::grid_for_each_neighbor(
                self.positions,
                grid,
                position,
                cutoff * cutoff,
                &mut emit,
            ),
            Some(CellGrid::Periodic(grid)) => {
                let periodic = self.periodic.as_ref().ok_or(SpatialError::InvalidCell)?;
                grid.for_each_neighbor(
                    self.positions,
                    position,
                    cutoff * cutoff,
                    periodic,
                    &mut emit,
                )
            }
            None => Ok(()),
        }
    }
}

/// Validates a non-negative finite cutoff.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_cutoff(cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

/// Validates a query cutoff against the cutoff used to construct the grid.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_query_cutoff(cutoff: f32, build_cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 && cutoff <= build_cutoff {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

/// Validates atom indices against a coordinate-array length.
///
/// Runtime is `O(N)` and no allocation is performed.
fn validate_indices(atoms: &[u32], position_count: usize) -> Result<(), SpatialError> {
    for &atom in atoms {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if index >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom));
        }
    }

    Ok(())
}

#[cfg(test)]
#[path = "cell_tests.rs"]
mod tests;
