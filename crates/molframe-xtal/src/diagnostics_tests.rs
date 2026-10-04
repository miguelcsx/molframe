use super::*;
use molframe_core::Diagnostic;

#[test]
fn map_and_reflection_failures_name_their_kind() {
    assert_eq!(Diagnostic::from(MrcError::Truncated).code(), Code::E7101);
    assert_eq!(
        Diagnostic::from(MrcError::UnsupportedMode(3)).code(),
        Code::E4105
    );
    assert_eq!(
        Diagnostic::from(MrcError::MemoryLimit {
            required: 2,
            limit: 1
        })
        .code(),
        Code::E7001
    );
    assert_eq!(
        Diagnostic::from(MapStatisticsError::MaskLength).code(),
        Code::E5102
    );
    assert_eq!(
        Diagnostic::from(ReflectionError::TruncatedMtz).code(),
        Code::E7101
    );
    assert_eq!(
        Diagnostic::from(ReflectionError::MillerIndices).code(),
        Code::E2004
    );
    assert!(
        Diagnostic::from(ReflectionError::InvalidMtz)
            .code()
            .is_registered()
    );
}
