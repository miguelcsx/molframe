//! Neighbour-joining phylogeny reconstruction.

mod algorithm;

pub use algorithm::{
    NegativeBranches, NjError, NjOptions, neighbor_joining, neighbor_joining_with,
    validate_distance_matrix,
};
