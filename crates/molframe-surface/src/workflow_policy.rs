//! Governed molecular-surface construction and intrinsic analysis.

use crate::{
    IndexedSurfaceMesh, SasaError, SurfaceCurvature, SurfaceDistances, SurfaceGeometryError,
    SurfaceGridOptions, edge_geodesic_distances, solvent_excluded_surface_with_options,
    surface_curvatures,
};
use molframe_core::contract::{AlgorithmId, Analysis, AnalysisPolicy, Coverage, ParameterValue};

/// Indexed mesh, curvature and optional source geodesics from one workflow.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceWorkflowResult {
    /// Validated shared-vertex surface.
    pub mesh: IndexedSurfaceMesh,
    /// One quality-labelled curvature estimate per vertex.
    pub curvature: Box<[SurfaceCurvature]>,
    /// Edge-geodesic distances when a source was requested.
    pub geodesics: Option<SurfaceDistances>,
}

/// Named scientific and numerical choices for a governed surface workflow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceWorkflowOptions<'a> {
    /// Stable atom-selection expression or identifier.
    pub selection: &'a str,
    /// Stable versioned radii-set identifier.
    pub radii_set: &'a str,
    /// Solvent probe radius in ångström.
    pub probe: f32,
    /// Surface-grid resolution and hard allocation ceiling.
    pub grid: SurfaceGridOptions,
    /// Optional mesh vertex from which to compute edge geodesics.
    pub source_vertex: Option<u32>,
}

/// Why a governed surface workflow was refused.
#[derive(Debug, thiserror::Error)]
pub enum SurfaceWorkflowError {
    /// Surface construction input or allocation was invalid.
    #[error(transparent)]
    Surface(#[from] SasaError),
    /// Intrinsic mesh analysis was invalid.
    #[error(transparent)]
    Geometry(#[from] SurfaceGeometryError),
    /// The public coverage counters cannot represent the atom count.
    #[error("surface input count exceeds u32 coverage limits")]
    CoverageOverflow,
}

/// Constructs an SES mesh and governed intrinsic descriptors.
///
/// `selection` and `radii_set` are stable caller-provided identifiers; the
/// numerical arrays remain the computation inputs.
///
/// # Errors
///
/// Returns a surface or mesh-geometry validation error.
pub fn governed_surface_geometry(
    positions: &[[f32; 3]],
    radii: &[f32],
    options: SurfaceWorkflowOptions<'_>,
    policy: &AnalysisPolicy,
) -> Result<Analysis<SurfaceWorkflowResult>, SurfaceWorkflowError> {
    let Some(probe) = ParameterValue::finite_float(f64::from(options.probe)) else {
        return Err(SasaError::InvalidProbe.into());
    };
    let Some(resolution) = ParameterValue::finite_float(f64::from(options.grid.resolution)) else {
        return Err(SasaError::InvalidGridOptions.into());
    };
    let max_grid_cells =
        i64::try_from(options.grid.max_cells).map_err(|_| SasaError::InvalidGridOptions)?;
    let surface =
        solvent_excluded_surface_with_options(positions, radii, options.probe, options.grid)?;
    let mesh = surface.indexed_mesh();
    let curvature = surface_curvatures(&mesh)?;
    let geodesics = match options.source_vertex {
        Some(source) => Some(edge_geodesic_distances(&mesh, source)?),
        None => None,
    };
    let intended =
        u32::try_from(positions.len()).map_err(|_| SurfaceWorkflowError::CoverageOverflow)?;
    let mut analysis = Analysis::complete(
        SurfaceWorkflowResult {
            mesh,
            curvature,
            geodesics,
        },
        Coverage::complete(intended),
        policy,
    );
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new("ses-cotangent-curvature", "1"))
        .with_parameter("selection", ParameterValue::Text(options.selection.into()))
        .with_parameter("radii_set", ParameterValue::Text(options.radii_set.into()))
        .with_parameter("probe", probe)
        .with_parameter("resolution", resolution)
        .with_parameter("max_grid_cells", ParameterValue::Integer(max_grid_cells))
        .with_parameter(
            "curvature_estimator",
            ParameterValue::Text("cotangent-angle-deficit".into()),
        );
    if let Some(source) = options.source_vertex {
        analysis.provenance = analysis
            .provenance
            .with_parameter("source_vertex", ParameterValue::Integer(i64::from(source)));
    }
    Ok(analysis)
}
