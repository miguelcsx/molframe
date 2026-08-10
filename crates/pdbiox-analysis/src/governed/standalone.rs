//! Governed entry points for scientific domains without a `Structure` topology.

use super::common::{descriptor, float, integer};
use crate::{
    FragmentMappingError, FragmentMatch, FragmentReference, PolymerError, PolymerStatistics,
    Pucker, map_fragments, polymer_statistics, sugar_pucker,
};
use pdbiox_core::contract::{
    Analysis, AnalysisPolicy, Coverage, MissingPolicy, ParameterValue, Provenance, SourceRef,
    Status,
};

/// A pure-domain kernel or governance failure.
#[derive(Debug, thiserror::Error)]
pub enum StandaloneAnalysisError<E> {
    /// The scientific kernel refused its inputs.
    #[error(transparent)]
    Kernel(E),
    /// Coverage cannot be represented in public counters.
    #[error("analysis input count exceeds u32 coverage limits")]
    CoverageOverflow,
    /// Policy requires failure for absent domain entries.
    #[error("analysis policy rejects {missing} missing inputs")]
    MissingData {
        /// Missing input entries.
        missing: u32,
    },
    /// This crate does not understand a newer missing-data policy.
    #[error("unsupported missing-data policy")]
    UnsupportedMissingPolicy,
}

fn coverage(intended: usize, used: usize) -> Option<Coverage> {
    let intended = u32::try_from(intended).ok()?;
    let used = u32::try_from(used).ok()?;
    Some(Coverage {
        intended,
        used,
        missing: intended - used,
        ambiguous: 0,
    })
}

fn governed<T, E>(
    value: T,
    coverage: Coverage,
    policy: &AnalysisPolicy,
    descriptor: &super::AnalysisDescriptor,
) -> Result<Analysis<T>, StandaloneAnalysisError<E>> {
    let status = if coverage.missing == 0 {
        Status::Complete
    } else {
        match policy.missing_atoms {
            MissingPolicy::Ignore | MissingPolicy::Report => Status::Partial,
            MissingPolicy::Indeterminate => Status::Indeterminate,
            MissingPolicy::Fail => {
                return Err(StandaloneAnalysisError::MissingData {
                    missing: coverage.missing,
                });
            }
            _ => return Err(StandaloneAnalysisError::UnsupportedMissingPolicy),
        }
    };
    Ok(Analysis {
        value,
        status,
        coverage,
        warnings: Vec::new(),
        assumptions: Vec::new(),
        provenance: descriptor.apply(Provenance::new(policy).with_source(SourceRef::Memory)),
    })
}

/// Maps an optional trace to an explicit fragment library under missing policy.
///
/// # Errors
///
/// Returns a library/input error, counter overflow, or missing-policy refusal.
pub fn governed_fragment_mapping(
    trace: &[Option<[f32; 3]>],
    library: &[FragmentReference],
    maximum_rmsd: f64,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<FragmentMatch>>, StandaloneAnalysisError<FragmentMappingError>> {
    let value =
        map_fragments(trace, library, maximum_rmsd).map_err(StandaloneAnalysisError::Kernel)?;
    let used = trace.iter().filter(|point| point.is_some()).count();
    let coverage = coverage(trace.len(), used).ok_or(StandaloneAnalysisError::CoverageOverflow)?;
    governed(
        value,
        coverage,
        policy,
        &descriptor("fragment-mapping")
            .with_parameter("library_size", integer(library.len()))
            .with_parameter("maximum_rmsd", float(maximum_rmsd)),
    )
}

/// Computes governed statistics for an explicitly ordered polymer path.
///
/// # Errors
///
/// Returns a polymer-domain error or coverage overflow.
pub fn governed_polymer_statistics(
    path: &[[f32; 3]],
    policy: &AnalysisPolicy,
) -> Result<Analysis<PolymerStatistics>, StandaloneAnalysisError<PolymerError>> {
    let value = polymer_statistics(path).map_err(StandaloneAnalysisError::Kernel)?;
    let coverage =
        coverage(path.len(), path.len()).ok_or(StandaloneAnalysisError::CoverageOverflow)?;
    governed(value, coverage, policy, &descriptor("polymer-statistics"))
}

/// Computes governed five-torsion sugar pseudorotation.
#[must_use]
pub fn governed_sugar_pucker(nu: [f64; 5], policy: &AnalysisPolicy) -> Analysis<Pucker> {
    Analysis {
        value: sugar_pucker(nu),
        status: Status::Complete,
        coverage: Coverage::complete(5),
        warnings: Vec::new(),
        assumptions: Vec::new(),
        provenance: descriptor("sugar-pucker")
            .with_parameter("torsions", ParameterValue::Text(format!("{nu:?}").into()))
            .apply(Provenance::new(policy).with_source(SourceRef::Memory)),
    }
}
