//! Contiguous fixed-radius cell list.

use crate::brute::{canonicalise, distance_squared, finite};
use crate::{CellGridOptions, NeighborPair, PeriodicBox, SpatialError};

#[path = "cell/geometry.rs"]
mod geometry;
use geometry::{finite_bounds, grid_geometry};
#[path = "cell/periodic.rs"]
mod periodic_grid;
use periodic_grid::PeriodicGrid;
#[path = "cell/unique.rs"]
mod unique;
pub(crate) use unique::{
    for_each_pairs_same_selection_unordered, pairs_same_selection, pairs_same_selection_unordered,
};

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

#[derive(Debug)]
struct Grid {
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    offsets: Vec<usize>,
    members: Vec<u32>,
    non_empty_cells: Vec<usize>,
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
    periodic: Option<&'a PeriodicBox>,
    grid: Option<CellGrid>,
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
        periodic: Option<&'a PeriodicBox>,
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
        periodic: Option<&'a PeriodicBox>,
        options: CellGridOptions,
    ) -> Result<Self, SpatialError> {
        validate_cutoff(cutoff)?;
        validate_indices(targets, positions.len())?;
        options.validate()?;

        let grid = match periodic {
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
        })
    }

    /// Finds pairs for the indexed targets.
    ///
    /// # Errors
    ///
    /// Returns an error when the requested cutoff exceeds the build cutoff or
    /// a query index is out of range.
    pub fn pairs(&self, query: &[u32], cutoff: f32) -> Result<Vec<NeighborPair>, SpatialError> {
        validate_query_cutoff(cutoff, self.cutoff)?;
        validate_indices(query, self.positions.len())?;

        match &self.grid {
            Some(CellGrid::Cartesian(grid)) => {
                grid_pairs(self.positions, grid, query, cutoff * cutoff)
            }
            Some(CellGrid::Periodic(grid)) => {
                let Some(periodic) = self.periodic else {
                    return Err(SpatialError::InvalidCell);
                };
                grid.pairs(self.positions, query, cutoff * cutoff, periodic)
            }
            None => Ok(Vec::new()),
        }
    }
}

/// Evaluates a cell-grid query and canonicalises the resulting pairs.
///
/// Each finite query atom visits at most 27 cells. Pairwise work inside those
/// cells depends on local spatial density.
fn grid_pairs(
    positions: &[[f32; 3]],
    grid: &Grid,
    query: &[u32],
    cutoff_squared: f32,
) -> Result<Vec<NeighborPair>, SpatialError> {
    let mut found = Vec::with_capacity(query.len());

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

        append_grid_pairs(positions, grid, atom, position, cutoff_squared, &mut found)?;
    }

    canonicalise(&mut found);
    Ok(found)
}

/// Searches the 27 neighbouring cells for one query atom.
///
/// Runtime is proportional to the members of the visited cells and no
/// temporary heap allocation is performed.
fn append_grid_pairs(
    positions: &[[f32; 3]],
    grid: &Grid,
    atom: u32,
    position: [f32; 3],
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) -> Result<(), SpatialError> {
    let centre = cell_of(grid, position)?;

    for delta in NEIGHBOUR_OFFSETS {
        let Some(cell) = neighbour_cell(grid.dims, centre, delta) else {
            continue;
        };

        let members = members_of_cell(grid, cell).ok_or(SpatialError::NumericRangeExceeded)?;
        append_cell_pairs(positions, atom, position, members, cutoff_squared, found)?;
    }
    Ok(())
}

/// Compares one query atom against all members of one cell.
///
/// Grid members are finite by construction, so only defensive bounds checks
/// remain in the hot loop.
fn append_cell_pairs(
    positions: &[[f32; 3]],
    atom: u32,
    position: [f32; 3],
    targets: &[u32],
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) -> Result<(), SpatialError> {
    for &target in targets {
        if atom == target {
            continue;
        }

        let Ok(index) = usize::try_from(target) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(target_position) = positions.get(index).copied() else {
            continue;
        };

        let squared = distance_squared(position, target_position, None);

        if squared <= cutoff_squared {
            found.push(NeighborPair::new(atom, target, squared));
        }
    }
    Ok(())
}

/// Builds the contiguous grid representation for finite target coordinates.
///
/// Building is `O(T + C)` for `T` targets and `C` cells and uses contiguous
/// count, offset and member arrays.
fn build_grid(
    positions: &[[f32; 3]],
    targets: &[u32],
    cutoff: f32,
    options: CellGridOptions,
) -> Result<Option<Grid>, SpatialError> {
    let Some(bounds) = finite_bounds(positions, targets) else {
        return Ok(None);
    };
    let (edge, dims, count) = if cutoff == 0.0 {
        (1.0, [1; 3], 1)
    } else {
        grid_geometry(bounds.min, bounds.max, cutoff, options)?
    };

    let counts = count_members(positions, targets, bounds.min, edge, dims, count)?;
    let offsets = prefix_offsets(&counts)?;
    let non_empty_cells = counts
        .iter()
        .enumerate()
        .filter_map(|(cell, &count)| (count != 0).then_some(cell))
        .collect();
    let members = fill_members(
        positions,
        targets,
        bounds.min,
        edge,
        dims,
        &offsets,
        bounds.finite_count,
    )?;

    Ok(Some(Grid {
        origin: bounds.min,
        edge,
        dims,
        offsets,
        members,
        non_empty_cells,
    }))
}

/// Computes one count per grid cell.
///
/// Runtime is `O(T)` with `O(C)` count storage.
fn count_members(
    positions: &[[f32; 3]],
    targets: &[u32],
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    count: usize,
) -> Result<Vec<usize>, SpatialError> {
    let mut counts = vec![0usize; count];

    for &target in targets {
        let Ok(index) = usize::try_from(target) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        if !finite(position) {
            continue;
        }

        let cell = linear(dims, raw_cell(origin, edge, dims, position)?);

        if let Some(slot) = counts.get_mut(cell) {
            *slot += 1;
        }
    }

    Ok(counts)
}

/// Converts per-cell counts into CSR-style prefix offsets.
///
/// Runtime is `O(C)` and output space is `O(C)`.
fn prefix_offsets(counts: &[usize]) -> Result<Vec<usize>, SpatialError> {
    let capacity = counts
        .len()
        .checked_add(1)
        .ok_or(SpatialError::NumericRangeExceeded)?;
    let mut offsets = Vec::with_capacity(capacity);
    let mut total = 0usize;

    offsets.push(total);

    for &count in counts {
        total = total
            .checked_add(count)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        offsets.push(total);
    }

    Ok(offsets)
}

/// Fills the contiguous CSR-style member array.
///
/// Runtime is `O(T + C)` and output space is exactly one `u32` per finite
/// target occurrence.
fn fill_members(
    positions: &[[f32; 3]],
    targets: &[u32],
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    offsets: &[usize],
    finite_count: usize,
) -> Result<Vec<u32>, SpatialError> {
    let cell_count = offsets
        .len()
        .checked_sub(1)
        .ok_or(SpatialError::NumericRangeExceeded)?;
    let mut cursors: Vec<usize> = offsets.iter().take(cell_count).copied().collect();
    let mut members = vec![0u32; finite_count];

    for &target in targets {
        let Ok(index) = usize::try_from(target) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        if !finite(position) {
            continue;
        }

        let cell = linear(dims, raw_cell(origin, edge, dims, position)?);

        let Some(cursor) = cursors.get_mut(cell) else {
            continue;
        };

        let Some(slot) = members.get_mut(*cursor) else {
            continue;
        };

        *slot = target;
        *cursor += 1;
    }

    Ok(members)
}

/// Converts a Cartesian point into clamped cell coordinates.
///
/// Runtime and auxiliary space are `O(1)`.
fn raw_cell(
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    point: [f32; 3],
) -> Result<[usize; 3], SpatialError> {
    let convert = |axis: usize| {
        let scaled = f64::from((point[axis] - origin[axis]) / edge);
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
fn cell_of(grid: &Grid, point: [f32; 3]) -> Result<[usize; 3], SpatialError> {
    raw_cell(grid.origin, grid.edge, grid.dims, point)
}

/// Converts three-dimensional cell coordinates to a linear cell index.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn linear(dims: [usize; 3], cell: [usize; 3]) -> usize {
    (cell[0] * dims[1] + cell[1]) * dims[2] + cell[2]
}

/// Returns a neighbouring linear cell index when it remains in bounds.
///
/// Runtime and auxiliary space are `O(1)`.
fn neighbour_cell(dims: [usize; 3], centre: [usize; 3], delta: [isize; 3]) -> Option<usize> {
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
fn members_of_cell(grid: &Grid, cell: usize) -> Option<&[u32]> {
    let &start = grid.offsets.get(cell)?;
    let next = cell.checked_add(1)?;
    let &end = grid.offsets.get(next)?;

    grid.members.get(start..end)
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
