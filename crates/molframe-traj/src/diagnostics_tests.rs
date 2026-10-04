use super::*;

#[test]
fn stream_and_analysis_failures_name_their_kind() {
    let mismatch = TrajectoryError::AtomCountMismatch {
        expected: 3,
        found: 4,
    };
    assert_eq!(Diagnostic::from(&mismatch).code(), Code::E5102);
    assert!(Diagnostic::from(mismatch).code().is_registered());
    assert_eq!(
        Diagnostic::from(TrajectoryError::Cancelled).code(),
        Code::E1904
    );
    assert_eq!(
        Diagnostic::from(TrajectoryError::MemoryLimit {
            required: 10,
            limit: 5
        })
        .code(),
        Code::E7001
    );
    assert_eq!(
        Diagnostic::from(TrajectoryError::SourceIo {
            format: "xtc",
            kind: std::io::ErrorKind::UnexpectedEof
        })
        .code(),
        Code::E7101
    );
}

#[test]
fn wrapped_failures_report_the_wrapped_code() {
    let build = TrajectoryError::Build(TrajectoryBuildError::DimensionOverflow);
    assert_eq!(Diagnostic::from(build).code(), Code::E1903);
    let spatial = TrajectoryError::Spatial(molframe_spatial::SpatialError::InvalidCell);
    assert_eq!(Diagnostic::from(spatial).code(), Code::E5004);
}

#[test]
fn ensemble_failures_name_their_kind() {
    assert_eq!(
        Diagnostic::from(EnsembleGeometryError::MemoryLimit).code(),
        Code::E7001
    );
    assert_eq!(
        Diagnostic::from(EnsembleGeometryError::Empty).code(),
        Code::E5103
    );
    assert_eq!(
        Diagnostic::from(KMeansError::DidNotConverge).code(),
        Code::E5104
    );
    assert_eq!(
        Diagnostic::from(EnsembleStatisticsError::PartialBlock).code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(MsdError::AtomOutOfBounds(9)).code(),
        Code::E6009
    );
    assert_eq!(
        Diagnostic::from(PathSimilarityError::InvalidPaths).code(),
        Code::E5102
    );
    assert_eq!(
        Diagnostic::from(FrameViewError::DimensionMismatch).code(),
        Code::E5102
    );
}

#[test]
fn a_wrapped_ensemble_failure_reports_the_wrapped_code() {
    let wrapped = GovernedEnsembleError::Similarity(EnsembleSimilarityError::SingularCovariance);
    assert_eq!(Diagnostic::from(wrapped).code(), Code::E5104);
    assert_eq!(
        Diagnostic::from(GovernedEnsembleError::CoverageOverflow).code(),
        Code::E1903
    );
}
