//! Solvent-accessible and molecular surfaces.
//!
//! These take positions and per-atom radii and return areas. They do not know
//! what an element is or which radius set produced the numbers they are given;
//! choosing the radii is a chemistry question answered a layer above, so that
//! the same surface code serves whatever convention a caller trusts.

#![forbid(unsafe_code)]

#[path = "sasa.rs"]
mod accessible_area;
#[path = "buried.rs"]
mod burial;
pub mod cavity;
mod cavity_geometry;
mod components;
pub mod depth;
pub mod geodesic;
mod grid_options;
mod io;
pub mod lee_richards;
#[path = "curvature.rs"]
mod local_shape;
mod neighbourhood;
mod numeric;
#[path = "sphere.rs"]
mod sampling;
pub mod ses;
#[path = "mesh.rs"]
mod triangulation;
#[path = "governed.rs"]
mod workflow_policy;

pub use accessible_area::{
    AtomContactArea, ExcludedSurfacePoint, SasaError, SurfacePoint, atom_contact_areas,
    shrake_rupley, surface_points, surface_points_at_density, surface_points_excluding_pairs,
};
pub use burial::{
    BuriedSurface, BuriedSurfaceError, MoleculeRole, buried_solvent_excluded_surface,
    buried_solvent_excluded_surface_with_options, buried_surface,
};
pub use cavity::{Cavity, cavities, cavities_with_options};
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
pub use sampling::fibonacci_sphere;
pub use ses::{
    SolventExcludedSurface, SurfaceTriangle, solvent_excluded_surface,
    solvent_excluded_surface_with_options,
};
pub use triangulation::{IndexedSurfaceMesh, MeshReport, SurfaceFace};
pub use workflow_policy::{
    SurfaceWorkflowError, SurfaceWorkflowOptions, SurfaceWorkflowResult, governed_surface_geometry,
};
