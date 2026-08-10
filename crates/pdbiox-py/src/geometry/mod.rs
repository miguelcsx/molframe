//! Vectorized geometry and rigid-transform bindings.

mod arrays;
mod intrinsic;
mod moments;
mod primitives;
mod rigid;

pub(crate) use arrays::{coordinates, distance_matrix, rmsd, superpose};
pub(crate) use intrinsic::{
    PyBackboneCoordinates, PyBackboneFrame, PyBackboneTorsions, PyCircularSummary, PyHelixGeometry,
    backbone_frames, circular_summary, helix_geometry, path_torsions, rotation_mean, torus_summary,
};
pub(crate) use moments::{
    PyAxes, PyEigenOptions, PyPlane, asphericity, best_fit_plane, centre_of_mass, centroid,
    gyration_axes, inertia_tensor, plane_deviation, principal_axes, radius_of_gyration,
};
pub(crate) use primitives::{
    angle, cross, degrees, dihedral, displacement, distance, distance_matrix_between,
    distance_squared, dot, norm, normalise, rmsf,
};
pub(crate) use rigid::{PyRigid, PySuperposition};
