use super::*;
use molframe_core::Diagnostic;

#[test]
fn each_failure_names_the_code_of_its_kind() {
    assert_eq!(
        Diagnostic::from(FluctuationError::NoFrames).code(),
        Code::E5103
    );
    assert_eq!(
        Diagnostic::from(FluctuationError::RaggedFrames).code(),
        Code::E5102
    );
    assert_eq!(
        Diagnostic::from(RotationError::DidNotConverge).code(),
        Code::E5104
    );
    assert_eq!(
        Diagnostic::from(PeriodicError::DimensionMismatch).code(),
        Code::E5102
    );
}
