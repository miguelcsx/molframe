//! Geometric kernels.
//!
//! Distances, angles, torsions, moments and superposition. Everything here is
//! deterministic and allocation-free in the measuring path, and everything that
//! sums over a set does so in double precision regardless of the precision the
//! positions were stored in.
//!
//! Nothing here knows what a residue is. These take positions and return
//! numbers; deciding *which* positions is a question for the layers above.

#![forbid(unsafe_code)]

/// Discrete geometry of ordered backbone traces.
pub mod backbone;
mod batch_measure;
mod diagnostics;
/// Deterministic symmetric eigendecomposition.
pub mod eigen;
/// Per-atom positional fluctuation across frames.
pub mod fluctuation;
/// Pairwise distance matrices.
pub mod matrix;
pub mod measure;
pub mod moments;
mod numeric;
/// Geometry along ordered coordinate paths.
pub mod path;
pub mod periodic_angle;
pub mod planar;
pub mod polymer;
/// Validated three-dimensional rotations.
pub mod rotation;
mod simd;
mod simd_measure;
/// Three-vector arithmetic in single precision.
pub mod single;
/// Rigid least-squares superposition.
pub mod superpose;
pub mod transform;

pub use backbone::{
    BackboneFrame, HelixGeometry, backbone_frames, helix_geometry, helix_geometry_with_options,
};
pub use batch_measure::{BatchGeometryError, angles_into, distances_into, torsions_into};
pub use eigen::{Decomposition, EigenError, EigenOptions, symmetric, symmetric_with_options};
pub use fluctuation::{FluctuationError, RmsfAccumulator, rmsf};
pub use matrix::{
    DistanceMatrix, MatrixError, condensed_distances_f32_with_context,
    condensed_distances_with_context, distance_block_rows, distance_matrix,
    distance_matrix_between, distance_matrix_f32_with_context, distance_matrix_into,
    distance_matrix_with_context, distance_rows_f32_into, distance_rows_into,
};
pub use measure::{
    angle, cross, degrees, dihedral, displacement, distance, distance_squared, dot, norm, normalise,
};
pub use moments::{
    SpreadFrame, asphericity, asphericity_with_options, centre_of_mass, centroid, gyration_axes,
    gyration_axes_with_options, inertia_tensor, principal_axes, principal_axes_with_options,
    radius_of_gyration, shape_parameter, shape_parameter_with_options, spread_frame,
};
pub use path::path_torsions;
pub use periodic_angle::{
    CircularSummary, PeriodicAngle, PeriodicError, TorusMetric, circular_summary, torus_summary,
};
pub use planar::{
    Plane, best_fit_plane, best_fit_plane_with_options, plane_deviation,
    plane_deviation_with_options,
};
pub use polymer::{BackboneResidue, BackboneTorsions, backbone_torsions};
pub use rotation::{
    Rotation3, RotationError, RotationMeanOptions, RotationOptions, rotation_mean,
    rotation_mean_with_options,
};
pub use single::{cross_f32, dot_f32, normalise_f32, sub_f32};
pub use superpose::{
    SuperposeError, SuperposeOptions, Superposition, rmsd, rmsd_after_fit, rmsd_flat, superpose,
    superpose_weighted, superpose_weighted_with_options, superpose_with_options,
};
pub use transform::Rigid;
