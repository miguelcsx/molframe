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

#[test]
fn mesh_and_burial_failures_name_their_kind() {
    assert_eq!(
        Diagnostic::from(SurfaceGeometryError::NonManifoldMesh).code(),
        Code::E5104
    );
    assert_eq!(
        Diagnostic::from(SurfaceGeometryError::VertexOutOfBounds).code(),
        Code::E5102
    );
    assert_eq!(
        Diagnostic::from(SurfaceComponentError::InvalidFilter).code(),
        Code::E5101
    );
    assert_eq!(Diagnostic::from(AtomDepthError).code(), Code::E5101);
    let wrapped = BuriedSurfaceError::Surface(SasaError::InvalidProbe);
    assert_eq!(Diagnostic::from(wrapped).code(), Code::E5101);
    let mismatch = BuriedSurfaceError::LengthMismatch {
        positions: 3,
        radii: 2,
        roles: 3,
    };
    assert_eq!(Diagnostic::from(mismatch).code(), Code::E5102);
}
