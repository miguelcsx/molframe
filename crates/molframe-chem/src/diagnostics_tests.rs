use super::*;
use molframe_core::Diagnostic;

#[test]
fn a_malformed_record_and_a_bad_pattern_name_their_kind() {
    assert_eq!(Diagnostic::from(MolError::Malformed).code(), Code::E1201);
    assert_eq!(Diagnostic::from(MolError::BondEndpoint).code(), Code::E3006);
    assert_eq!(
        Diagnostic::from(Mol2Error::CountMismatch).code(),
        Code::E5102
    );
    let pattern = SmartsError {
        position: 4,
        message: "unclosed bracket".into(),
    };
    let diagnostic = Diagnostic::from(&pattern);
    assert_eq!(diagnostic.code(), Code::E1301);
    assert_eq!(
        Diagnostic::from(SmartsDataError::Connectivity).code(),
        Code::E4003
    );
    assert_eq!(
        Diagnostic::from(AutomorphismLimit { limit: 10 }).code(),
        Code::E1901
    );
}
