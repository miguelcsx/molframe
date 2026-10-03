use super::*;
use molframe_core::Diagnostic;

#[test]
fn alignment_failures_name_their_kind_and_keep_the_message() {
    let none = Diagnostic::from(AlignError::NoAlignmentPath);
    assert_eq!(none.code(), Code::E5105);
    assert!(none.code().is_registered());
    assert_eq!(
        Diagnostic::from(AlignError::NumericOverflow).code(),
        Code::E1903
    );
    assert_eq!(
        Diagnostic::from(RegionAlignError::Alignment(AlignError::NoAlignmentPath)).code(),
        Code::E5105
    );
}

#[test]
fn sequence_content_failures_share_the_data_validity_code() {
    let unknown = SequenceError::UnknownSymbol {
        position: 3,
        symbol: b'!',
    };
    let diagnostic = Diagnostic::from(unknown);
    assert_eq!(diagnostic.code(), Code::E2101);
    assert!(
        Diagnostic::from(MatrixError::UnknownProfile)
            .code()
            .is_registered()
    );
    assert_eq!(
        Diagnostic::from(NewickError::UnexpectedEnd).code(),
        Code::E2101
    );
}

#[test]
fn a_memory_ceiling_is_the_resource_code() {
    let error = MsaError::MemoryLimit {
        required: 10,
        limit: 5,
    };
    assert_eq!(Diagnostic::from(error).code(), Code::E7001);
}
