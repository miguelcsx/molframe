//! Contiguous fixed-radius cell list.

use crate::brute::{canonicalise, distance_squared, finite};
use crate::{NeighborPair, PeriodicBox, SpatialError};

const MAX_CELLS: usize = 1_000_000;

#[derive(Debug)]
struct Grid {
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    offsets: Vec<usize>,
    members: Vec<u32>,
}

/// A fixed-radius grid with one contiguous member array.
///
/// Building is O(n); a query visits at most twenty-seven cells when the index
/// is used at or below its construction cutoff.
#[derive(Debug)]
pub struct CellList<'a> {
    positions: &'a [[f32; 3]],
    targets: Vec<u32>,
    cutoff: f32,
    periodic: Option<&'a PeriodicBox>,
    grid: Option<Grid>,
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
        if !cutoff.is_finite() || cutoff < 0.0 {
            return Err(SpatialError::InvalidCutoff);
        }
        for target in targets {
            if *target as usize >= positions.len() {
                return Err(SpatialError::AtomOutOfBounds(*target));
            }
        }
        let grid = if periodic.is_some() || cutoff == 0.0 {
            None
        } else {
            build_grid(positions, targets, cutoff)
        };
        Ok(Self {
            positions,
            targets: targets.to_vec(),
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
        if !cutoff.is_finite() || cutoff < 0.0 || cutoff > self.cutoff {
            return Err(SpatialError::InvalidCutoff);
        }
        for atom in query {
            if *atom as usize >= self.positions.len() {
                return Err(SpatialError::AtomOutOfBounds(*atom));
            }
        }
        let Some(grid) = &self.grid else {
            return Ok(crate::brute::pairs(
                self.positions,
                query,
                &self.targets,
                cutoff * cutoff,
                self.periodic,
            ));
        };
        let cutoff_squared = cutoff * cutoff;
        let mut found = Vec::new();
        for atom in query {
            let Some(position) = self.positions.get(*atom as usize).copied() else {
                continue;
            };
            if !finite(position) {
                continue;
            }
            let centre = cell_of(grid, position);
            for di in -1..=1 {
                for dj in -1..=1 {
                    for dk in -1..=1 {
                        let Some(cell) = neighbour_cell(grid.dims, centre, [di, dj, dk]) else {
                            continue;
                        };
                        let range = grid.offsets[cell]..grid.offsets[cell + 1];
                        for target in &grid.members[range] {
                            if atom == target {
                                continue;
                            }
                            let Some(target_position) =
                                self.positions.get(*target as usize).copied()
                            else {
                                continue;
                            };
                            let squared = distance_squared(position, target_position, None);
                            if squared <= cutoff_squared {
                                found.push(NeighborPair::new(*atom, *target, squared));
                            }
                        }
                    }
                }
            }
        }
        canonicalise(&mut found);
        Ok(found)
    }
}

fn build_grid(positions: &[[f32; 3]], targets: &[u32], cutoff: f32) -> Option<Grid> {
    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    let mut finite_count = 0usize;
    for target in targets {
        let position = *positions.get(*target as usize)?;
        if !finite(position) {
            continue;
        }
        finite_count += 1;
        for axis in 0..3 {
            min[axis] = min[axis].min(position[axis]);
            max[axis] = max[axis].max(position[axis]);
        }
    }
    if finite_count == 0 {
        return None;
    }

    let mut edge = cutoff;
    let mut dims = dimensions(min, max, edge);
    while cell_count(dims) > MAX_CELLS {
        edge *= 2.0;
        dims = dimensions(min, max, edge);
    }
    let count = cell_count(dims);
    let mut counts = vec![0usize; count];
    for target in targets {
        let Some(position) = positions.get(*target as usize).copied() else {
            continue;
        };
        if finite(position) {
            counts[linear(dims, raw_cell(min, edge, dims, position))] += 1;
        }
    }
    let mut offsets = Vec::with_capacity(count + 1);
    offsets.push(0);
    for count in counts {
        let next = offsets.last().copied().map_or(0, |offset| offset + count);
        offsets.push(next);
    }
    let mut cursors = offsets[..count].to_vec();
    let mut members = vec![0u32; finite_count];
    for target in targets {
        let Some(position) = positions.get(*target as usize).copied() else {
            continue;
        };
        if !finite(position) {
            continue;
        }
        let cell = linear(dims, raw_cell(min, edge, dims, position));
        if let Some(cursor) = cursors.get_mut(cell)
            && let Some(slot) = members.get_mut(*cursor)
        {
            *slot = *target;
            *cursor += 1;
        }
    }
    Some(Grid {
        origin: min,
        edge,
        dims,
        offsets,
        members,
    })
}

fn dimensions(min: [f32; 3], max: [f32; 3], edge: f32) -> [usize; 3] {
    std::array::from_fn(|axis| (((max[axis] - min[axis]) / edge).floor() as usize + 1).max(1))
}

fn cell_count(dims: [usize; 3]) -> usize {
    dims[0].saturating_mul(dims[1]).saturating_mul(dims[2])
}

fn raw_cell(origin: [f32; 3], edge: f32, dims: [usize; 3], point: [f32; 3]) -> [usize; 3] {
    std::array::from_fn(|axis| {
        (((point[axis] - origin[axis]) / edge).floor() as usize).min(dims[axis] - 1)
    })
}

fn cell_of(grid: &Grid, point: [f32; 3]) -> [usize; 3] {
    raw_cell(grid.origin, grid.edge, grid.dims, point)
}

fn linear(dims: [usize; 3], cell: [usize; 3]) -> usize {
    (cell[0] * dims[1] + cell[1]) * dims[2] + cell[2]
}

fn neighbour_cell(dims: [usize; 3], centre: [usize; 3], delta: [i32; 3]) -> Option<usize> {
    let mut cell = [0usize; 3];
    for axis in 0..3 {
        let signed = isize::try_from(delta[axis]).ok()?;
        let coordinate = centre[axis].checked_add_signed(signed)?;
        if coordinate >= dims[axis] {
            return None;
        }
        cell[axis] = coordinate;
    }
    Some(linear(dims, cell))
}

#[cfg(test)]
#[path = "cell_tests.rs"]
mod tests;
