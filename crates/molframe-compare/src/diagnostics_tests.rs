use super::*;

#[test]
fn comparison_failures_name_the_kind_of_input_that_was_wrong() {
    let length = CompareError::LengthMismatch {
        model: 3,
        reference: 4,
    };
    assert_eq!(Diagnostic::from(&length).code(), Code::E5102);
    assert_eq!(
        Diagnostic::from(CompareError::InvalidDistanceCutoff).code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(CompareError::NoComparablePairs).code(),
        Code::E5003
    );
    assert!(Diagnostic::from(length).code().is_registered());
}

#[test]
fn a_governed_failure_reports_the_code_of_the_failure_it_wraps() {
    let wrapped = GovernedCompareError::Compare(CompareError::InvalidScoreInput);
    assert_eq!(Diagnostic::from(&wrapped).code(), Code::E5101);
    let carried = GovernedCompareError::Diagnostic(Diagnostic::new(Code::E4003));
    assert_eq!(Diagnostic::from(carried).code(), Code::E4003);
    assert_eq!(
        Diagnostic::from(GovernedCompareError::UnsupportedPolicy("assembly")).code(),
        Code::E6103
    );
}

#[test]
fn a_mapping_failure_keeps_its_own_diagnostic() {
    let inner = Diagnostic::new(Code::E4003);
    let error = MappedCompareError::Mapping(inner.clone());
    assert_eq!(Diagnostic::from(error).code(), inner.code());
    let scored = MappedCompareError::Compare(CompareError::InvalidMapping);
    assert_eq!(Diagnostic::from(scored).code(), Code::E5102);
}
