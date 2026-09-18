use super::PyCifValue;

#[test]
fn cif_value_projection_preserves_sentinels_and_numeric_kinds() {
    let unknown = PyCifValue::from(molframe::CifValue::Unknown);
    assert_eq!(unknown.kind(), "unknown");
    assert!(!unknown.is_recorded());

    let integer = PyCifValue::from(molframe::CifValue::Integer(42));
    assert_eq!(integer.kind(), "integer");
    assert_eq!(integer.integer(), Some(42));
    assert_eq!(integer.number(), Some(42.0));
    assert_eq!(integer.as_identifier().as_deref(), Some("42"));
}
