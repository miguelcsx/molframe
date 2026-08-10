//! Per-group ensemble dispersion and block convergence summaries.

use crate::numeric::f64_from_usize;

/// Coordinate dispersion for one explicit atom group.
#[derive(Clone, Debug, PartialEq)]
pub struct GroupVariance {
    /// Atom indices used for the group, in caller order.
    pub atoms: Vec<usize>,
    /// Mean population variance along x, y, and z over group atoms.
    pub variance_by_axis: [f64; 3],
    /// Square root of the sum of the three mean axis variances.
    pub rms_fluctuation: f64,
}

/// One complete or caller-permitted partial convergence block.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConvergenceBlock {
    /// Inclusive first observation index.
    pub start: usize,
    /// Exclusive last observation index.
    pub end: usize,
    /// Mean within this block.
    pub block_mean: f64,
    /// Mean from observation zero through this block.
    pub cumulative_mean: f64,
}

/// Policy for a final block smaller than the requested size.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RemainderPolicy {
    /// Include the final partial block and report its exact bounds.
    Include,
    /// Reject a series whose length is not a multiple of the block size.
    Reject,
}

/// Why ensemble statistics could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum EnsembleStatisticsError {
    /// Frames must form a non-empty finite rectangular coordinate ensemble.
    #[error("frames must form a non-empty finite rectangular coordinate ensemble")]
    InvalidFrames,
    /// Groups must be non-empty and contain only valid atom indices.
    #[error("groups must be non-empty and contain valid atom indices")]
    InvalidGroups,
    /// Scalar observations and block size must be valid and finite.
    #[error("values must be non-empty and finite and block size non-zero")]
    InvalidBlocks,
    /// Remainder was rejected by explicit policy.
    #[error("observation count is not a multiple of the requested block size")]
    PartialBlock,
}

/// Computes coordinate variance for explicit atom groups over pre-aligned frames.
///
/// Frames are not fitted implicitly. Callers choose and apply any alignment
/// transform before requesting dispersion, preserving the meaning of variance.
///
/// # Errors
///
/// Returns [`EnsembleStatisticsError`] for invalid frames or groups.
pub fn group_coordinate_variance(
    frames: &[Vec<[f32; 3]>],
    groups: &[Vec<usize>],
) -> Result<Vec<GroupVariance>, EnsembleStatisticsError> {
    let atom_count = validate_frames(frames)?;
    if groups.is_empty()
        || groups
            .iter()
            .any(|group| group.is_empty() || group.iter().any(|atom| *atom >= atom_count))
    {
        return Err(EnsembleStatisticsError::InvalidGroups);
    }
    groups
        .iter()
        .map(|group| group_variance(frames, group))
        .collect()
}

/// Summarizes scalar convergence using an explicit block and remainder policy.
///
/// # Errors
///
/// Returns [`EnsembleStatisticsError`] for invalid values, block size, or remainder.
pub fn block_convergence(
    values: &[f64],
    block_size: usize,
    remainder: RemainderPolicy,
) -> Result<Vec<ConvergenceBlock>, EnsembleStatisticsError> {
    if values.is_empty() || block_size == 0 || values.iter().any(|value| !value.is_finite()) {
        return Err(EnsembleStatisticsError::InvalidBlocks);
    }
    if remainder == RemainderPolicy::Reject && !values.len().is_multiple_of(block_size) {
        return Err(EnsembleStatisticsError::PartialBlock);
    }
    let mut cumulative_sum = 0.0;
    values
        .chunks(block_size)
        .enumerate()
        .map(|(block, chunk)| {
            let start = block * block_size;
            let end = start + chunk.len();
            let block_sum: f64 = chunk.iter().sum();
            cumulative_sum += block_sum;
            let block_count =
                f64_from_usize(chunk.len()).ok_or(EnsembleStatisticsError::InvalidBlocks)?;
            let cumulative_count =
                f64_from_usize(end).ok_or(EnsembleStatisticsError::InvalidBlocks)?;
            Ok(ConvergenceBlock {
                start,
                end,
                block_mean: block_sum / block_count,
                cumulative_mean: cumulative_sum / cumulative_count,
            })
        })
        .collect()
}

fn validate_frames(frames: &[Vec<[f32; 3]>]) -> Result<usize, EnsembleStatisticsError> {
    let Some(first) = frames.first() else {
        return Err(EnsembleStatisticsError::InvalidFrames);
    };
    if first.is_empty()
        || frames.iter().any(|frame| frame.len() != first.len())
        || frames
            .iter()
            .flatten()
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(EnsembleStatisticsError::InvalidFrames);
    }
    Ok(first.len())
}

fn group_variance(
    frames: &[Vec<[f32; 3]>],
    group: &[usize],
) -> Result<GroupVariance, EnsembleStatisticsError> {
    let observations = frames
        .len()
        .checked_mul(group.len())
        .and_then(f64_from_usize)
        .ok_or(EnsembleStatisticsError::InvalidFrames)?;
    let means: [f64; 3] = std::array::from_fn(|axis| {
        frames
            .iter()
            .flat_map(|frame| group.iter().map(|atom| f64::from(frame[*atom][axis])))
            .sum::<f64>()
            / observations
    });
    let variance_by_axis = std::array::from_fn(|axis| {
        frames
            .iter()
            .flat_map(|frame| group.iter().map(|atom| f64::from(frame[*atom][axis])))
            .map(|value| (value - means[axis]).powi(2))
            .sum::<f64>()
            / observations
    });
    Ok(GroupVariance {
        atoms: group.to_vec(),
        rms_fluctuation: variance_by_axis.into_iter().sum::<f64>().sqrt(),
        variance_by_axis,
    })
}

#[cfg(test)]
#[path = "ensemble_statistics_tests.rs"]
mod tests;
