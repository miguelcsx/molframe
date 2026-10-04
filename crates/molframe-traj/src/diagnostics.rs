//! Stable diagnostic codes for every error the trajectory layer returns.

use crate::{
    EnsembleGeometryError, EnsembleSimilarityError, EnsembleStatisticsError, FrameViewError,
    GovernedEnsembleError, KMeansError, MsdError, PathSimilarityError, TrajectoryBuildError,
    TrajectoryError,
};
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

diagnostic_from!(FrameViewError, |error| match error {
    FrameViewError::DimensionOverflow => Code::E1903,
    FrameViewError::DimensionMismatch => Code::E5102,
});

diagnostic_from!(EnsembleGeometryError, |error| match error {
    EnsembleGeometryError::Empty => Code::E5103,
    EnsembleGeometryError::DimensionMismatch => Code::E5102,
    EnsembleGeometryError::InvalidParameter => Code::E5101,
    EnsembleGeometryError::MemoryLimit => Code::E7001,
    EnsembleGeometryError::DegenerateFit | EnsembleGeometryError::DidNotConverge => Code::E5104,
});

diagnostic_from!(KMeansError, |error| match error {
    KMeansError::InvalidObservations | KMeansError::InvalidOptions => Code::E5101,
    KMeansError::EmptyCluster(_) | KMeansError::DidNotConverge => Code::E5104,
});

diagnostic_from!(EnsembleStatisticsError, |_error| Code::E5101);

diagnostic_from!(EnsembleSimilarityError, |error| match error {
    EnsembleSimilarityError::InvalidEnsembles => Code::E5102,
    EnsembleSimilarityError::InvalidOptions | EnsembleSimilarityError::InvalidLabels => {
        Code::E5101
    }
    EnsembleSimilarityError::MemoryLimit => Code::E7001,
    EnsembleSimilarityError::SingularCovariance => Code::E5104,
});

diagnostic_from!(MsdError, |error| match error {
    MsdError::Empty => Code::E5103,
    MsdError::DimensionMismatch => Code::E5102,
    MsdError::AtomOutOfBounds(_) => Code::E6009,
    MsdError::NonFinite | MsdError::InvalidLag => Code::E5101,
    MsdError::ObservationLimit => Code::E1903,
});

diagnostic_from!(PathSimilarityError, |error| match error {
    PathSimilarityError::InvalidPaths => Code::E5102,
    PathSimilarityError::MemoryLimit => Code::E7001,
    PathSimilarityError::DegenerateFit => Code::E5104,
});

impl From<&GovernedEnsembleError> for Diagnostic {
    /// A wrapped failure names the code of the failure it wraps.
    fn from(error: &GovernedEnsembleError) -> Self {
        match error {
            GovernedEnsembleError::Geometry(inner) => Self::from(inner),
            GovernedEnsembleError::Similarity(inner) => Self::from(inner),
            GovernedEnsembleError::KMeans(inner) => Self::from(inner),
            GovernedEnsembleError::Statistics(inner) => Self::from(inner),
            GovernedEnsembleError::MeanSquaredDisplacement(inner) => Self::from(inner),
            GovernedEnsembleError::CoverageOverflow => {
                Self::new(Code::E1903).with_message(error.to_string())
            }
        }
    }
}

impl From<GovernedEnsembleError> for Diagnostic {
    fn from(error: GovernedEnsembleError) -> Self {
        Self::from(&error)
    }
}

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
