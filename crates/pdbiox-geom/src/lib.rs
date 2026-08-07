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
// Vectors have components called x, y, z and quaternions have one called w.
// Renaming them to satisfy a length rule would make this code harder to check
// against the mathematics it implements, which is the only way to check it.
#![allow(clippy::many_single_char_names, clippy::similar_names)]

pub mod eigen;
pub mod matrix;
pub mod measure;
pub mod moments;
pub mod path;
pub mod polymer;
pub mod superpose;
pub mod transform;

pub use matrix::{DistanceMatrix, distance_matrix, distance_matrix_between};
pub use measure::{
    angle, cross, degrees, dihedral, displacement, distance, distance_squared, dot, norm, normalise,
};
pub use moments::{
    asphericity, centre_of_mass, centroid, gyration_axes, inertia_tensor, principal_axes,
    radius_of_gyration,
};
pub use path::path_torsions;
pub use polymer::{BackboneResidue, BackboneTorsions, backbone_torsions};
pub use superpose::{SuperposeError, Superposition, rmsd, superpose};
pub use transform::Rigid;
