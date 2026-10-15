//! Solvent-accessible and molecular surfaces.
//!
//! These take positions and per-atom radii and return areas. They do not know
//! what an element is or which radius set produced the numbers they are given;
//! choosing the radii is a chemistry question answered a layer above, so that
//! the same surface code serves whatever convention a caller trusts.

#![forbid(unsafe_code)]

mod diagnostics;
mod exposure;
mod grid_options;
mod io;
#[cfg(not(target_arch = "wasm32"))]
pub mod native_solver;
mod numeric;
mod solver_output;
mod workflow_policy;
mod workspace;
pub use exposure::lee_richards;
use exposure::{
    accessible_area, burial, collect, neighbourhood, sampling, sampling_plan, stream,
    structure_sasa,
};
mod shape;
pub use shape::cavity;
pub use shape::depth;
pub use shape::geodesic;
pub use shape::pocket;
pub use shape::ses;
use shape::{cavity_geometry, components, local_shape, mesh_diagnostics, triangulation};

pub use accessible_area::{
    AtomContactArea, ExcludedSurfacePoint, SasaError, SurfacePoint, atom_contact_areas,
    shrake_rupley, surface_points, surface_points_at_density, surface_points_excluding_pairs,
};
pub use burial::{
    BuriedSurface, BuriedSurfaceError, MoleculeRole, MoleculeRoleError,
    buried_solvent_excluded_surface, buried_solvent_excluded_surface_with_options, buried_surface,
};
pub use cavity::{Cavity, cavities, cavities_with_options};
pub use collect::collect_shrake_rupley;
pub use components::{
    SurfaceComponent, SurfaceComponentError, SurfaceComponentFilter, filter_surface_components,
    surface_components,
};
pub use depth::{AtomDepthError, AtomDepthOptions, atom_depths};
pub use geodesic::{
    SurfaceDistances, SurfaceGeometryError, edge_geodesic_distances, surface_patch,
};
pub use grid_options::SurfaceGridOptions;
pub use io::write_obj;
pub use lee_richards::lee_richards;
pub use local_shape::{CurvatureQuality, SurfaceCurvature, surface_curvatures};
pub use pocket::{GridPocket, PocketScanError, PocketScanOptions, PointProximity, scan_pockets};
pub use sampling::fibonacci_sphere;
pub use sampling_plan::SasaSampler;
pub use ses::{
    SolventExcludedSurface, SurfaceTriangle, solvent_excluded_surface,
    solvent_excluded_surface_with_options,
};
pub use solver_output::{
    MolecularSurfaceMesh, SurfaceOutputError, SurfaceOutputLimits, read_msms_surface,
    read_surf_surface,
};
pub use stream::visit_shrake_rupley;
pub use structure_sasa::{
    DEFAULT_POINTS, DEFAULT_PROBE, DEFAULT_SLICES, SasaMethod, SasaOptions, SasaReport,
    SasaWordError, measure_atoms, structure_atom_areas, structure_sasa_report,
};
pub use triangulation::{IndexedSurfaceMesh, MeshReport, SurfaceFace};
pub use workflow_policy::{
    SurfaceWorkflowError, SurfaceWorkflowOptions, SurfaceWorkflowResult, governed_surface_geometry,
};
