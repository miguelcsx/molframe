use super::*;

#[test]
fn wrapped_errors_name_the_code_of_what_they_wrap() {
    assert_eq!(
        Diagnostic::from(PlanarityError::InvalidTolerance).code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(BFactorKernelError::Validation(BFactorError::NoObservations)).code(),
        Code::E5103
    );
    assert_eq!(
        Diagnostic::from(CompletenessError::ExplicitNamespace).code(),
        Code::E6103
    );
    assert_eq!(
        Diagnostic::from(RealSpaceCorrelationError::ZeroVariance).code(),
        Code::E5104
    );
}
