//! Bounded voxel-grid geometry construction.

use crate::{SasaError, numeric::f64_to_usize};

pub(crate) fn geometry(
    positions: &[[f32; 3]],
    expanded: &[f64],
    step: f64,
    max_cells: usize,
) -> Result<([f64; 3], [usize; 3]), SasaError> {
    let (low, high) = spatial_bounds(positions, expanded);
    let mut origin = [0.0; 3];
    let mut dims = [0usize; 3];
    let mut total = 1usize;
    for axis in 0..3 {
        origin[axis] = low[axis] - step;
        let span = (high[axis] + step) - origin[axis];
        let count = f64_to_usize((span / step).ceil())
            .checked_add(1)
            .ok_or(SasaError::GridDimensionsOverflow)?
            .max(3);
        dims[axis] = count;
        total = total
            .checked_mul(count)
            .ok_or(SasaError::GridDimensionsOverflow)?;
        if total > max_cells {
            return Err(SasaError::GridTooLarge { cells: total });
        }
    }
    Ok((origin, dims))
}

fn spatial_bounds(positions: &[[f32; 3]], expanded: &[f64]) -> ([f64; 3], [f64; 3]) {
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for (position, &radius) in positions.iter().zip(expanded) {
        for axis in 0..3 {
            let centre = f64::from(position[axis]);
            low[axis] = low[axis].min(centre - radius);
            high[axis] = high[axis].max(centre + radius);
        }
    }
    (low, high)
}
