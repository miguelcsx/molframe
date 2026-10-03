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
