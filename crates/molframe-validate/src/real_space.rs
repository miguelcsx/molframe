//! Deterministic real-space Pearson correlation for density maps.

use molframe_xtal::{DensityMap, MapBoundary};

use crate::numeric::usize_to_f64;

/// Summary of a real-space map correlation calculation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealSpaceCorrelation {
    /// Pearson product-moment correlation coefficient.
    pub coefficient: f64,
    /// Number of paired density samples included.
    pub sample_count: usize,
    /// Mean density in the observed map over included samples.
    pub observed_mean: f64,
    /// Mean density in the calculated map over included samples.
    pub calculated_mean: f64,
}

/// Failure to define a real-space correlation coefficient.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum RealSpaceCorrelationError {
    /// Grids or a supplied mask do not describe the same voxel set.
    #[error("real-space correlation inputs have incompatible grids or mask")]
    IncompatibleGrid,
    /// Fewer than two paired finite samples were available.
    #[error("real-space correlation requires at least two paired samples")]
    InsufficientSamples,
    /// At least one paired density was non-finite.
    #[error("real-space correlation input contains a non-finite density")]
    NonFiniteDensity,
    /// One of the sampled fields has no variance.
    #[error("real-space correlation is undefined for a constant density field")]
    ZeroVariance,
}

/// Correlates all corresponding voxels of two identically placed maps.
///
/// # Errors
///
/// Returns an error for incompatible map grids, non-finite densities, fewer
/// than two samples, or a constant field.
pub fn real_space_map_correlation(
    observed: &DensityMap,
    calculated: &DensityMap,
) -> Result<RealSpaceCorrelation, RealSpaceCorrelationError> {
    ensure_compatible(observed, calculated)?;
    correlate_pairs(
        observed
            .values
            .iter()
            .copied()
            .zip(calculated.values.iter().copied()),
    )
}

/// Correlates the selected voxels of two identically placed maps.
///
/// # Errors
///
/// Returns an error for incompatible grids or mask length, non-finite
/// densities, fewer than two selected samples, or a constant selected field.
pub fn masked_real_space_correlation(
    observed: &DensityMap,
    calculated: &DensityMap,
    mask: &[bool],
) -> Result<RealSpaceCorrelation, RealSpaceCorrelationError> {
    ensure_compatible(observed, calculated)?;
    if mask.len() != observed.values.len() {
        return Err(RealSpaceCorrelationError::IncompatibleGrid);
    }
    correlate_pairs(
        observed
            .values
            .iter()
            .copied()
            .zip(calculated.values.iter().copied())
            .zip(mask)
            .filter_map(|(pair, include)| include.then_some(pair)),
    )
}

/// Samples two maps at Cartesian coordinates and correlates the paired values.
///
/// Coordinates missing from either non-periodic map are omitted and reflected
/// in `sample_count`.
///
/// # Errors
///
/// Returns an error for fewer than two shared samples, non-finite densities,
/// or a constant sampled field.
pub fn sampled_real_space_correlation(
    observed: &DensityMap,
    calculated: &DensityMap,
    positions: &[[f64; 3]],
    boundary: MapBoundary,
) -> Result<RealSpaceCorrelation, RealSpaceCorrelationError> {
    // Both maps build their cell transform once here rather than once per
    // sampled coordinate, which is where nearly all the time went.
    // A degenerate cell previously made every sample return nothing, which
    // reached the caller as an empty sample set; that outcome is preserved.
    let (Some(observed), Some(calculated)) = (observed.sampler(), calculated.sampler()) else {
        return Err(RealSpaceCorrelationError::InsufficientSamples);
    };

    correlate_pairs(positions.iter().filter_map(|position| {
        observed
            .sample_cartesian(*position, boundary)
            .zip(calculated.sample_cartesian(*position, boundary))
    }))
}

fn ensure_compatible(
    observed: &DensityMap,
    calculated: &DensityMap,
) -> Result<(), RealSpaceCorrelationError> {
    if observed.dimensions != calculated.dimensions
        || observed.starts != calculated.starts
        || observed.sampling != calculated.sampling
        || observed.cell != calculated.cell
        || !same_f64_triplet(observed.origin, calculated.origin)
        || observed.values.len() != calculated.values.len()
    {
        return Err(RealSpaceCorrelationError::IncompatibleGrid);
    }
    Ok(())
}

fn same_f64_triplet(left: [f64; 3], right: [f64; 3]) -> bool {
    left.into_iter()
        .zip(right)
        .all(|(left, right)| left.to_bits() == right.to_bits())
}

fn correlate_pairs<I>(pairs: I) -> Result<RealSpaceCorrelation, RealSpaceCorrelationError>
where
    I: IntoIterator<Item = (f32, f32)>,
{
    let mut count = 0_usize;
    let mut observed_mean = 0.0;
    let mut calculated_mean = 0.0;
    let mut observed_sum_squares = 0.0;
    let mut calculated_sum_squares = 0.0;
    let mut co_moment = 0.0;
    for (observed, calculated) in pairs {
        if !observed.is_finite() || !calculated.is_finite() {
            return Err(RealSpaceCorrelationError::NonFiniteDensity);
        }
        count += 1;
        let count_f64 = usize_to_f64(count);
        let observed_delta = f64::from(observed) - observed_mean;
        let calculated_delta = f64::from(calculated) - calculated_mean;
        observed_mean += observed_delta / count_f64;
        calculated_mean += calculated_delta / count_f64;
        observed_sum_squares += observed_delta * (f64::from(observed) - observed_mean);
        calculated_sum_squares += calculated_delta * (f64::from(calculated) - calculated_mean);
        co_moment += observed_delta * (f64::from(calculated) - calculated_mean);
    }
    if count < 2 {
        return Err(RealSpaceCorrelationError::InsufficientSamples);
    }
    let denominator = (observed_sum_squares * calculated_sum_squares).sqrt();
    if denominator <= f64::EPSILON {
        return Err(RealSpaceCorrelationError::ZeroVariance);
    }
    Ok(RealSpaceCorrelation {
        coefficient: co_moment / denominator,
        sample_count: count,
        observed_mean,
        calculated_mean,
    })
}

#[cfg(test)]
#[path = "real_space_tests.rs"]
mod tests;
