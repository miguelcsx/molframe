//! Vectorized geometry and rigid-transform bindings.

mod arrays;
mod distance_types;
mod intrinsic;
pub(crate) mod intrinsic_geometry;
mod moments;
mod primitives;
mod rigid;
mod types;

pub(crate) use arrays::{
    borrowed_coordinates, coordinates, distance_matrix, rmsd, rmsd_flat, superpose,
    superpose_with_options,
};
pub(crate) use distance_types::PyDistanceMatrix;
pub(crate) use intrinsic::{
    PyBackboneCoordinates, PyBackboneFrame, PyBackboneTorsions, PyCircularSummary, PyHelixGeometry,
    backbone_frames, circular_summary, helix_geometry, helix_geometry_with_options, path_torsions,
    rotation_mean, rotation_mean_with_options, torus_summary,
};
pub(crate) use moments::{
    PyAxes, PyEigenOptions, PyPlane, asphericity, asphericity_with_options, best_fit_plane,
    best_fit_plane_with_options, centre_of_mass, centroid, gyration_axes,
    gyration_axes_with_options, inertia_tensor, plane_deviation, plane_deviation_with_options,
    principal_axes, principal_axes_with_options, radius_of_gyration,
};
pub(crate) use primitives::{
    angle, cross, degrees, dihedral, displacement, distance, distance_matrix_between,
    distance_squared, dot, norm, normalise, rmsf,
};
pub(crate) use rigid::{PyRigid, PySuperposition};
pub(crate) use types::{
    PyRotationMeanOptions, PyRotationOptions, PySuperposeOptions, backbone_torsions, register,
    symmetric, symmetric_with_options,
};
