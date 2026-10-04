use super::*;
use crate::policy_execution::GovernedAnalysisError;

#[test]
fn a_kernel_failure_inside_the_governed_wrapper_keeps_its_own_code() {
    let wrapped = GovernedAnalysisError::Kernel(HseError::MissingRoleAnnotation);
    assert_eq!(
        governed_diagnostic(&wrapped, |inner| Diagnostic::from(inner)).code(),
        Code::E6103
    );
    let missing: GovernedAnalysisError<HseError> = GovernedAnalysisError::MissingData {
        missing: 2,
        ambiguous: 0,
    };
    assert_eq!(
        governed_diagnostic(&missing, |inner| Diagnostic::from(inner)).code(),
        Code::E5103
    );
}

#[test]
fn each_kind_of_failure_names_its_kind_of_code() {
    assert_eq!(
        Diagnostic::from(PiStackingError::InvalidOptions).code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(HydrogenBondError::MissingCell).code(),
        Code::E5004
    );
    assert_eq!(Diagnostic::from(PolymerError::TooShort).code(), Code::E5103);
    assert_eq!(
        Diagnostic::from(NativeError::AtomCountMismatch {
            reference: 1,
            target: 2
        })
        .code(),
        Code::E5102
    );
    assert_eq!(
        Diagnostic::from(NucleicTorsionError::MissingRoles).code(),
        Code::E6103
    );
}
