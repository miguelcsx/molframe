//! Stable diagnostic codes for every error the trajectory layer returns.

use crate::{TrajectoryBuildError, TrajectoryError};
use molframe_core::{Code, Diagnostic, diagnostic_from};

diagnostic_from!(TrajectoryBuildError, |error| match error {
    TrajectoryBuildError::DimensionOverflow => Code::E1903,
    TrajectoryBuildError::AtomCountMismatch { .. } => Code::E5102,
});

impl From<&TrajectoryError> for Diagnostic {
    /// Wrapped failures name the code of the failure they wrap; the rest name
    /// the kind of thing that went wrong with the stream or the analysis.
    fn from(error: &TrajectoryError) -> Self {
        let code = match error {
            TrajectoryError::Build(inner) => return Self::from(inner),
            TrajectoryError::Spatial(inner) => return Self::from(inner),
            TrajectoryError::RandomAccessUnavailable | TrajectoryError::AnalysisNotParallel => {
                Code::E6103
            }
            TrajectoryError::AtomCountMismatch { .. } => Code::E5102,
            TrajectoryError::SelectionOutOfRange { .. } => Code::E6009,
            TrajectoryError::DegenerateFit => Code::E5104,
            TrajectoryError::InvalidWorkerCount => Code::E6102,
            TrajectoryError::WorkerPanicked => Code::E9001,
            TrajectoryError::CoordinateGenerationExhausted
            | TrajectoryError::UnrepresentableCoordinate
            | TrajectoryError::IdentityOverflow => Code::E1903,
            TrajectoryError::MissingFrame { .. }
            | TrajectoryError::SourceIo { .. }
            | TrajectoryError::InvalidSource { .. } => Code::E7101,
            TrajectoryError::MissingCell => Code::E5004,
            TrajectoryError::OverlappingGroups { .. } => Code::E5101,
            TrajectoryError::MemoryLimit { .. } => Code::E7001,
            TrajectoryError::Cancelled => Code::E1904,
        };
        Self::new(code).with_message(error.to_string())
    }
}

impl From<TrajectoryError> for Diagnostic {
    fn from(error: TrajectoryError) -> Self {
        Self::from(&error)
    }
}

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
