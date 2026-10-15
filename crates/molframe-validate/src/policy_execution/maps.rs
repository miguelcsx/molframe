//! Governed real-space map correlation entry points.

use crate::{
    RealSpaceCorrelation, RealSpaceCorrelationError, masked_real_space_correlation,
    real_space_map_correlation, sampled_real_space_correlation,
};
use molframe_core::contract::{
    AlgorithmId, Analysis, AnalysisPolicy, Coverage, MissingPolicyError, PeriodicPolicy,
    Provenance, Quality, SourceRef, resolve_missing,
};
use molframe_xtal::{DensityMap, MapBoundary};

/// Real-space kernel, policy or coverage failure.
#[derive(Debug, thiserror::Error)]
pub enum GovernedMapError {
    /// The correlation kernel refused its maps or samples.
    #[error(transparent)]
    Correlation(#[from] RealSpaceCorrelationError),
    /// Public coverage counters cannot represent the sample count.
    #[error("real-space sample count exceeds u32 coverage limits")]
    CoverageOverflow,
    /// Missing samples are forbidden by policy.
    #[error("real-space correlation omitted {0} requested samples")]
    MissingSamples(u32),
    /// The installed policy has a variant this crate does not understand.
    #[error("unsupported analysis policy field: {0}")]
    UnsupportedPolicy(&'static str),
}

fn result(
    value: RealSpaceCorrelation,
    intended: usize,
    policy: &AnalysisPolicy,
    algorithm: &'static str,
) -> Result<Analysis<RealSpaceCorrelation>, GovernedMapError> {
    let intended = u32::try_from(intended).map_err(|_| GovernedMapError::CoverageOverflow)?;
    let used = u32::try_from(value.sample_count).map_err(|_| GovernedMapError::CoverageOverflow)?;
    let missing = intended
        .checked_sub(used)
        .ok_or(GovernedMapError::CoverageOverflow)?;
    let coverage = Coverage {
        intended,
        used,
        missing,
        ambiguous: 0,
    };
    let (outcome, quality) =
        resolve_missing(value, Quality::Complete, coverage, policy.missing_atoms).map_err(
            |error| match error {
                MissingPolicyError::Fail { missing, .. } => {
                    GovernedMapError::MissingSamples(missing)
                }
                _ => GovernedMapError::UnsupportedPolicy("missing_atoms"),
            },
        )?;
    Ok(Analysis::from_parts(
        outcome,
        quality,
        coverage,
        Vec::new(),
        Vec::new(),
        Provenance::new(policy)
            .with_source(SourceRef::Memory)
            .with_algorithm(AlgorithmId::new(algorithm, "1")),
    ))
}

/// Correlates every voxel of two compatible maps with governed provenance.
///
/// # Errors
///
/// Returns correlation or coverage errors.
pub fn governed_real_space_map_correlation(
    observed: &DensityMap,
    calculated: &DensityMap,
    policy: &AnalysisPolicy,
) -> Result<Analysis<RealSpaceCorrelation>, GovernedMapError> {
    let value = real_space_map_correlation(observed, calculated)?;
    result(
        value,
        observed.values.len(),
        policy,
        "real-space-correlation",
    )
}

/// Correlates an explicit voxel mask with governed provenance.
///
/// # Errors
///
/// Returns correlation, missing-policy or coverage errors.
pub fn governed_masked_real_space_correlation(
    observed: &DensityMap,
    calculated: &DensityMap,
    mask: &[bool],
    policy: &AnalysisPolicy,
) -> Result<Analysis<RealSpaceCorrelation>, GovernedMapError> {
    let value = masked_real_space_correlation(observed, calculated, mask)?;
    result(
        value,
        mask.iter().filter(|include| **include).count(),
        policy,
        "masked-real-space-correlation",
    )
}

/// Samples two maps at Cartesian positions under policy-driven boundaries.
///
/// # Errors
///
/// Returns correlation, policy, missing-policy or coverage errors.
pub fn governed_sampled_real_space_correlation(
    observed: &DensityMap,
    calculated: &DensityMap,
    positions: &[[f64; 3]],
    policy: &AnalysisPolicy,
) -> Result<Analysis<RealSpaceCorrelation>, GovernedMapError> {
    let boundary = match policy.periodic {
        PeriodicPolicy::None => MapBoundary::Missing,
        PeriodicPolicy::Pbc | PeriodicPolicy::MinimumImage => MapBoundary::Periodic,
        _ => return Err(GovernedMapError::UnsupportedPolicy("periodic")),
    };
    let value = sampled_real_space_correlation(observed, calculated, positions, boundary)?;
    result(
        value,
        positions.len(),
        policy,
        "sampled-real-space-correlation",
    )
}

#[cfg(test)]
#[path = "maps_tests.rs"]
mod tests;
