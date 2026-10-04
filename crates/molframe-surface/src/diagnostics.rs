//! Stable diagnostic codes for the surface kernels.

use crate::{
    AtomDepthError, BuriedSurfaceError, MoleculeRoleError, SasaError, SurfaceComponentError,
    SurfaceGeometryError,
};
use molframe_core::{Code, Diagnostic, diagnostic_from};

diagnostic_from!(SasaError, |error| match error {
    SasaError::LengthMismatch { .. } => Code::E5102,
    SasaError::InvalidProbe
    | SasaError::InvalidRadius
    | SasaError::InvalidDensity
    | SasaError::NoPoints
    | SasaError::InvalidGridOptions => Code::E5101,
    SasaError::GridTooLarge { .. }
    | SasaError::GridDimensionsOverflow
    | SasaError::AllocationFailed { .. } => Code::E1903,
    SasaError::WorkspaceTooLarge { .. } => Code::E7001,
    SasaError::Spatial(inner) => return Diagnostic::from(inner),
    SasaError::WorkerPanicked => Code::E9001,
});

diagnostic_from!(BuriedSurfaceError, |error| match error {
    BuriedSurfaceError::LengthMismatch { .. } => Code::E5102,
    BuriedSurfaceError::Surface(inner) => return Diagnostic::from(inner),
});

diagnostic_from!(SurfaceGeometryError, |error| match error {
    SurfaceGeometryError::VertexOutOfBounds => Code::E5102,
    SurfaceGeometryError::NonManifoldMesh => Code::E5104,
    SurfaceGeometryError::InvalidRadius => Code::E5101,
});

diagnostic_from!(SurfaceComponentError, |error| match error {
    SurfaceComponentError::InvalidFilter => Code::E5101,
    SurfaceComponentError::MeshTooLarge => Code::E1903,
});

diagnostic_from!(AtomDepthError, |_error| Code::E5101);

diagnostic_from!(MoleculeRoleError, |_error| Code::E5101);

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
