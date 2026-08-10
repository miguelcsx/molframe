//! Periodic image bounds and deterministic lattice traversal.

use crate::{KdPeriodicOptions, PeriodicBox, SpatialError};

pub(super) fn periodic_image_limits(
    periodic: &PeriodicBox,
    cutoff: f32,
) -> Result<[i64; 3], SpatialError> {
    let reach = periodic.fractional_cutoff_bounds(cutoff);
    let convert = |value| crate::numeric::ceil_i64(value + 1.0).ok_or(SpatialError::InvalidCell);
    Ok([convert(reach[0])?, convert(reach[1])?, convert(reach[2])?])
}

pub(super) fn validate_image_budget(
    limits: [i64; 3],
    options: KdPeriodicOptions,
) -> Result<(), SpatialError> {
    let required = limits.into_iter().try_fold(1usize, |count, limit| {
        let value = usize::try_from(limit).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let width = value
            .checked_mul(2)
            .and_then(|width| width.checked_add(1))
            .ok_or(SpatialError::NumericRangeExceeded)?;
        count
            .checked_mul(width)
            .ok_or(SpatialError::NumericRangeExceeded)
    })?;
    if required > options.maximum_image_count {
        Err(SpatialError::PeriodicImageLimitExceeded {
            required,
            maximum: options.maximum_image_count,
        })
    } else {
        Ok(())
    }
}

pub(super) fn for_each_shift(
    limits: [i64; 3],
    mut visit: impl FnMut([i64; 3]) -> Result<(), SpatialError>,
) -> Result<(), SpatialError> {
    for first in -limits[0]..=limits[0] {
        for second in -limits[1]..=limits[1] {
            for third in -limits[2]..=limits[2] {
                visit([first, second, third])?;
            }
        }
    }
    Ok(())
}
