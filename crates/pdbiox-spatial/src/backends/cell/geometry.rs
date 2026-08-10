//! Cell-grid bounds and memory-bounded geometry.

use crate::brute::finite;
use crate::numeric::floor_usize;
use crate::{CellGridOptions, SpatialError};

#[derive(Clone, Copy)]
pub(super) struct Bounds {
    pub(super) min: [f32; 3],
    pub(super) max: [f32; 3],
    pub(super) finite_count: usize,
}

pub(super) fn finite_bounds(positions: &[[f32; 3]], targets: &[u32]) -> Option<Bounds> {
    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    let mut finite_count = 0usize;

    for &target in targets {
        let index = usize::try_from(target).ok()?;
        let position = positions.get(index).copied()?;
        if !finite(position) {
            continue;
        }
        finite_count += 1;
        update_bounds(&mut min, &mut max, position);
    }

    (finite_count != 0).then_some(Bounds {
        min,
        max,
        finite_count,
    })
}

#[inline]
fn update_bounds(min: &mut [f32; 3], max: &mut [f32; 3], position: [f32; 3]) {
    for axis in 0..3 {
        min[axis] = min[axis].min(position[axis]);
        max[axis] = max[axis].max(position[axis]);
    }
}

pub(super) fn grid_geometry(
    min: [f32; 3],
    max: [f32; 3],
    cutoff: f32,
    options: CellGridOptions,
) -> Result<(f32, [usize; 3], usize), SpatialError> {
    let mut edge = cutoff;
    let mut dims = dimensions(min, max, edge)?;
    let mut count = cell_count(dims)?;

    while count > options.maximum_cell_count {
        edge *= options.edge_growth_factor;
        dims = dimensions(min, max, edge)?;
        count = cell_count(dims)?;
    }

    Ok((edge, dims, count))
}

fn dimensions(min: [f32; 3], max: [f32; 3], edge: f32) -> Result<[usize; 3], SpatialError> {
    Ok([
        dimension(min[0], max[0], edge)?,
        dimension(min[1], max[1], edge)?,
        dimension(min[2], max[2], edge)?,
    ])
}

fn dimension(min: f32, max: f32, edge: f32) -> Result<usize, SpatialError> {
    let cells =
        floor_usize(f64::from((max - min) / edge)).ok_or(SpatialError::NumericRangeExceeded)?;
    cells
        .checked_add(1)
        .map(|dimension| dimension.max(1))
        .ok_or(SpatialError::NumericRangeExceeded)
}

#[inline]
fn cell_count(dims: [usize; 3]) -> Result<usize, SpatialError> {
    dims[0]
        .checked_mul(dims[1])
        .and_then(|plane| plane.checked_mul(dims[2]))
        .ok_or(SpatialError::NumericRangeExceeded)
}
