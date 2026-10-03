use super::*;

#[test]
fn surface_failures_name_their_kind() {
    assert_eq!(
        Diagnostic::from(SasaError::InvalidProbe).code(),
        Code::E5101
    );
    assert_eq!(
        Diagnostic::from(SasaError::LengthMismatch {
            positions: 3,
            radii: 2
        })
        .code(),
        Code::E5102
    );
    assert_eq!(
        Diagnostic::from(SasaError::WorkspaceTooLarge { bytes: 9, limit: 1 }).code(),
        Code::E7001
    );
    let wrapped = SasaError::Spatial(molframe_spatial::SpatialError::InvalidCell);
    assert_eq!(Diagnostic::from(wrapped).code(), Code::E5004);
}
