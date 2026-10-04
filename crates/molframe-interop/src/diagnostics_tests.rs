use super::*;

#[test]
fn every_failure_names_the_code_of_its_kind() {
    assert_eq!(
        Diagnostic::from(GraphError::MissingCell).code(),
        Code::E5004
    );
    assert_eq!(
        Diagnostic::from(GraphError::ResourceLimit).code(),
        Code::E7001
    );
    assert_eq!(
        Diagnostic::from(GraphError::UnsupportedFeature { feature: "element" }).code(),
        Code::E6103
    );
    assert_eq!(
        Diagnostic::from(DatasetError::IndexOutOfBounds { index: 4 }).code(),
        Code::E6009
    );
    assert_eq!(
        Diagnostic::from(DatasetError::InvalidBatchSize).code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(DatasetError::DuplicateId { id: "1abc".into() }).code(),
        Code::E7101
    );
    assert_eq!(
        Diagnostic::from(DlpackError::RaggedCoordinates).code(),
        Code::E5102
    );
}
