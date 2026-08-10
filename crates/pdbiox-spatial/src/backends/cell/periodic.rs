//! Fractional-coordinate cell grid on a triclinic periodic torus.

use crate::brute::{canonicalise, distance_squared, finite};
use crate::{CellGridOptions, NeighborPair, PeriodicBox, SpatialError};

const NEIGHBOUR_OFFSETS: [[isize; 3]; 27] = super::NEIGHBOUR_OFFSETS;

#[derive(Debug)]
pub(super) struct PeriodicGrid {
    dims: [usize; 3],
    offsets: Vec<usize>,
    members: Vec<u32>,
}

impl PeriodicGrid {
    pub(super) fn build(
        positions: &[[f32; 3]],
        targets: &[u32],
        cutoff: f32,
        periodic: &PeriodicBox,
        options: CellGridOptions,
    ) -> Result<Self, SpatialError> {
        let dims = periodic_dimensions(periodic, cutoff, options)?;
        let cell_count = cell_count(dims)?;
        let mut counts = vec![0usize; cell_count];

        for &target in targets {
            let index = usize::try_from(target).map_err(|_| SpatialError::NumericRangeExceeded)?;
            let Some(position) = positions.get(index).copied() else {
                continue;
            };
            if finite(position) {
                counts[linear(dims, cell(periodic, position, dims)?)?] += 1;
            }
        }

        let offsets = prefix_offsets(&counts)?;
        let mut cursors: Vec<usize> = offsets.iter().take(cell_count).copied().collect();
        let finite_count = match offsets.last() {
            Some(count) => *count,
            None => 0,
        };
        let mut members = vec![0u32; finite_count];

        for &target in targets {
            let index = usize::try_from(target).map_err(|_| SpatialError::NumericRangeExceeded)?;
            let Some(position) = positions.get(index).copied() else {
                continue;
            };
            if !finite(position) {
                continue;
            }
            let index = linear(dims, cell(periodic, position, dims)?)?;
            let Some(cursor) = cursors.get_mut(index) else {
                continue;
            };
            if let Some(slot) = members.get_mut(*cursor) {
                *slot = target;
                *cursor += 1;
            }
        }

        Ok(Self {
            dims,
            offsets,
            members,
        })
    }

    pub(super) fn pairs(
        &self,
        positions: &[[f32; 3]],
        query: &[u32],
        cutoff_squared: f32,
        periodic: &PeriodicBox,
    ) -> Result<Vec<NeighborPair>, SpatialError> {
        let mut found = Vec::new();

        for &atom in query {
            let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
            let Some(position) = positions.get(index).copied() else {
                continue;
            };
            if !finite(position) {
                continue;
            }
            self.append_pairs(
                positions,
                atom,
                position,
                cutoff_squared,
                periodic,
                &mut found,
            )?;
        }

        canonicalise(&mut found);
        Ok(found)
    }

    fn append_pairs(
        &self,
        positions: &[[f32; 3]],
        atom: u32,
        position: [f32; 3],
        cutoff_squared: f32,
        periodic: &PeriodicBox,
        found: &mut Vec<NeighborPair>,
    ) -> Result<(), SpatialError> {
        let centre = cell(periodic, position, self.dims)?;
        let mut visited = [usize::MAX; NEIGHBOUR_OFFSETS.len()];
        let mut visited_count = 0usize;

        for delta in NEIGHBOUR_OFFSETS {
            let candidate = periodic_neighbour(self.dims, centre, delta)?;
            if visited[..visited_count].contains(&candidate) {
                continue;
            }
            visited[visited_count] = candidate;
            visited_count += 1;

            let members = self
                .members(candidate)
                .ok_or(SpatialError::NumericRangeExceeded)?;
            for &target in members {
                if atom == target {
                    continue;
                }
                let index =
                    usize::try_from(target).map_err(|_| SpatialError::NumericRangeExceeded)?;
                let Some(target_position) = positions.get(index).copied() else {
                    continue;
                };
                let squared = distance_squared(position, target_position, Some(periodic));
                if squared <= cutoff_squared {
                    found.push(NeighborPair::new(atom, target, squared));
                }
            }
        }
        Ok(())
    }

    fn members(&self, cell: usize) -> Option<&[u32]> {
        let &start = self.offsets.get(cell)?;
        let next = cell.checked_add(1)?;
        let &end = self.offsets.get(next)?;
        self.members.get(start..end)
    }
}

fn periodic_dimensions(
    periodic: &PeriodicBox,
    cutoff: f32,
    options: CellGridOptions,
) -> Result<[usize; 3], SpatialError> {
    if cutoff == 0.0 {
        return Ok([1; 3]);
    }
    let mut reach = periodic.fractional_cutoff_bounds(cutoff);
    let mut dims = dimensions(reach, options.maximum_cell_count)?;

    while cell_count(dims)? > options.maximum_cell_count {
        reach = reach.map(|value| value * f64::from(options.edge_growth_factor));
        dims = dimensions(reach, options.maximum_cell_count)?;
    }
    Ok(dims)
}

fn dimensions(reach: [f64; 3], maximum: usize) -> Result<[usize; 3], SpatialError> {
    Ok([
        dimension_for_reach(reach[0], maximum)?,
        dimension_for_reach(reach[1], maximum)?,
        dimension_for_reach(reach[2], maximum)?,
    ])
}

fn dimension_for_reach(reach: f64, maximum: usize) -> Result<usize, SpatialError> {
    if reach <= 0.0 {
        Ok(1)
    } else {
        let value =
            crate::numeric::floor_usize(reach.recip()).ok_or(SpatialError::NumericRangeExceeded)?;
        Ok(value.min(maximum).max(1))
    }
}

fn cell(
    periodic: &PeriodicBox,
    position: [f32; 3],
    dims: [usize; 3],
) -> Result<[usize; 3], SpatialError> {
    let fractional = periodic
        .fractional(position)
        .map(|value| value.rem_euclid(1.0));
    let coordinate = |axis: usize| {
        let dimension =
            crate::numeric::usize_f64(dims[axis]).ok_or(SpatialError::NumericRangeExceeded)?;
        let value = crate::numeric::floor_usize(fractional[axis] * dimension)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let maximum = dims[axis].checked_sub(1).ok_or(SpatialError::InvalidCell)?;
        Ok(value.min(maximum))
    };
    Ok([coordinate(0)?, coordinate(1)?, coordinate(2)?])
}

fn periodic_neighbour(
    dims: [usize; 3],
    centre: [usize; 3],
    delta: [isize; 3],
) -> Result<usize, SpatialError> {
    let coordinate = |axis| {
        let dimension =
            isize::try_from(dims[axis]).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let origin =
            isize::try_from(centre[axis]).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let shifted = origin
            .checked_add(delta[axis])
            .ok_or(SpatialError::NumericRangeExceeded)?;
        usize::try_from(shifted.rem_euclid(dimension))
            .map_err(|_| SpatialError::NumericRangeExceeded)
    };
    linear(dims, [coordinate(0)?, coordinate(1)?, coordinate(2)?])
}

#[inline]
fn linear(dims: [usize; 3], cell: [usize; 3]) -> Result<usize, SpatialError> {
    cell[0]
        .checked_mul(dims[1])
        .and_then(|value| value.checked_add(cell[1]))
        .and_then(|value| value.checked_mul(dims[2]))
        .and_then(|value| value.checked_add(cell[2]))
        .ok_or(SpatialError::NumericRangeExceeded)
}

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

fn cell_count(dims: [usize; 3]) -> Result<usize, SpatialError> {
    dims[0]
        .checked_mul(dims[1])
        .and_then(|value| value.checked_mul(dims[2]))
        .ok_or(SpatialError::NumericRangeExceeded)
}
