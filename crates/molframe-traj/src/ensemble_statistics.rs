//! Per-group ensemble dispersion and block convergence summaries.

use crate::frame_view::{FrameSource, FrameView};
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
    group_coordinate_variance_source(frames, groups)
}

/// Computes coordinate variance directly from borrowed contiguous frames.
///
/// The input coordinates remain borrowed; only the requested result records
/// are allocated.
///
/// # Errors
///
/// Returns [`EnsembleStatisticsError`] for invalid frames or groups.
pub fn group_coordinate_variance_view(
    frames: FrameView<'_>,
    groups: &[Vec<usize>],
) -> Result<Vec<GroupVariance>, EnsembleStatisticsError> {
    group_coordinate_variance_source(&frames, groups)
}

fn group_coordinate_variance_source<S: FrameSource + ?Sized>(
    frames: &S,
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
    // Once the complete input length is exactly representable, every block
    // boundary and block length is exactly representable too.  Validate that
    // invariant once instead of paying the checked conversion cost for every
    // output row.
    let total_count = f64_from_usize(values.len()).ok_or(EnsembleStatisticsError::InvalidBlocks)?;
    let full_block_count = (block_size <= values.len())
        .then(|| f64_from_usize(block_size))
        .flatten();
    let remainder_count = (!values.len().is_multiple_of(block_size))
        .then(|| f64_from_usize(values.len() % block_size))
        .flatten();
    let mut cumulative_sum = 0.0;
    let mut cumulative_count = 0.0;
    let mut blocks = Vec::with_capacity(values.len().div_ceil(block_size));
    for (block, chunk) in values.chunks(block_size).enumerate() {
        let start = block * block_size;
        let end = start + chunk.len();
        let block_sum: f64 = chunk.iter().sum();
        cumulative_sum += block_sum;
        let block_count = if chunk.len() == block_size {
            full_block_count.ok_or(EnsembleStatisticsError::InvalidBlocks)?
        } else {
            match remainder_count {
                Some(count) => count,
                None => total_count,
            }
        };
        cumulative_count += block_count;
        blocks.push(ConvergenceBlock {
            start,
            end,
            block_mean: block_sum / block_count,
            cumulative_mean: cumulative_sum / cumulative_count,
        });
    }
    Ok(blocks)
}

fn validate_frames<S: FrameSource + ?Sized>(frames: &S) -> Result<usize, EnsembleStatisticsError> {
    if frames.frame_count() == 0 {
        return Err(EnsembleStatisticsError::InvalidFrames);
    }
    let atom_count = frames.atom_count();
    if atom_count == 0 {
        return Err(EnsembleStatisticsError::InvalidFrames);
    }
    for index in 0..frames.frame_count() {
        let frame = frames
            .frame(index)
            .ok_or(EnsembleStatisticsError::InvalidFrames)?;
        if frame.len() != atom_count || frame.iter().flatten().any(|value| !value.is_finite()) {
            return Err(EnsembleStatisticsError::InvalidFrames);
        }
    }
    Ok(atom_count)
}

fn group_variance<S: FrameSource + ?Sized>(
    frames: &S,
    group: &[usize],
) -> Result<GroupVariance, EnsembleStatisticsError> {
    let observations = frames
        .frame_count()
        .checked_mul(group.len())
        .and_then(f64_from_usize)
        .ok_or(EnsembleStatisticsError::InvalidFrames)?;
    let mut sums = [0.0; 3];
    for frame_index in 0..frames.frame_count() {
        let frame = frames
            .frame(frame_index)
            .ok_or(EnsembleStatisticsError::InvalidFrames)?;
        for &atom in group {
            for axis in 0..3 {
                sums[axis] += f64::from(frame[atom][axis]);
            }
        }
    }
    let means = sums.map(|sum| sum / observations);
    let mut variance_by_axis = [0.0; 3];
    for frame_index in 0..frames.frame_count() {
        let frame = frames
            .frame(frame_index)
            .ok_or(EnsembleStatisticsError::InvalidFrames)?;
        for &atom in group {
            for axis in 0..3 {
                variance_by_axis[axis] += (f64::from(frame[atom][axis]) - means[axis]).powi(2);
            }
        }
    }
    variance_by_axis = variance_by_axis.map(|sum| sum / observations);
    Ok(GroupVariance {
        atoms: group.to_vec(),
        rms_fluctuation: variance_by_axis.into_iter().sum::<f64>().sqrt(),
        variance_by_axis,
    })
}

#[cfg(test)]
#[path = "ensemble_statistics_tests.rs"]
mod tests;
