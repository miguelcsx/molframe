//! Governed calibration-free nucleic frame parameters.

use super::common::{descriptor, float};
use super::standalone::StandaloneAnalysisError;
use crate::{
    BaseFrame, HelicalError, HelicalOptions, HelicalParameters, helical_parameters, helical_steps,
};
use molframe_core::contract::{Analysis, AnalysisPolicy, Coverage, Provenance, SourceRef, Status};

fn complete<T>(
    value: T,
    frame_count: usize,
    options: HelicalOptions,
    policy: &AnalysisPolicy,
    algorithm: &'static str,
) -> Result<Analysis<T>, StandaloneAnalysisError<HelicalError>> {
    let frame_count =
        u32::try_from(frame_count).map_err(|_| StandaloneAnalysisError::CoverageOverflow)?;
    let provenance = descriptor(algorithm)
        .with_parameter("frame_tolerance", float(options.frame_tolerance))
        .apply(Provenance::new(policy).with_source(SourceRef::Memory));
    Ok(Analysis {
        value,
        status: Status::Complete,
        coverage: Coverage::complete(frame_count),
        warnings: Vec::new(),
        assumptions: Vec::new(),
        provenance,
    })
}

/// Computes governed base-pair or step parameters for two explicit frames.
///
/// # Errors
///
/// Returns a frame-domain error or coverage overflow.
pub fn governed_helical_parameters(
    first: BaseFrame,
    second: BaseFrame,
    options: HelicalOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<HelicalParameters>, StandaloneAnalysisError<HelicalError>> {
    let value =
        helical_parameters(first, second, options).map_err(StandaloneAnalysisError::Kernel)?;
    complete(value, 2, options, policy, "nucleic-helical-parameters")
}

/// Computes governed parameters for consecutive explicit base frames.
///
/// # Errors
///
/// Returns a frame-domain error or coverage overflow.
pub fn governed_helical_steps(
    frames: &[BaseFrame],
    options: HelicalOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<HelicalParameters>>, StandaloneAnalysisError<HelicalError>> {
    let value = helical_steps(frames, options).map_err(StandaloneAnalysisError::Kernel)?;
    complete(
        value,
        frames.len(),
        options,
        policy,
        "nucleic-helical-steps",
    )
}
