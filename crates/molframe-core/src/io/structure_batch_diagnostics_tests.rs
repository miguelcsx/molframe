use super::*;

#[test]
fn a_wrapped_diagnostic_keeps_its_own_code_and_limits_name_theirs() {
    let inner = Diagnostic::new(Code::E7101).with_message("truncated");
    let wrapped = StructureBatchError::Diagnostic(inner.clone());
    assert_eq!(Diagnostic::from(&wrapped).code(), Code::E7101);
    assert_eq!(
        Diagnostic::from(StructureBatchError::DemandTooSmall {
            required: 9,
            available: 1
        })
        .code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(StructureBatchError::RecordExceedsBudget {
            required: 9,
            available: 1
        })
        .code(),
        Code::E7001
    );
    assert_eq!(
        Diagnostic::from(StructureBatchError::DictionaryFull).code(),
        Code::E1903
    );
}
