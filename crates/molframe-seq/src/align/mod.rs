//! Pairwise sequence alignment with affine gap costs.

mod api;
mod dynamic;
mod linear;
mod small;
mod types;

pub use api::{
    align_region, global, global_banded, global_matrix, global_score, local, local_matrix,
    semi_global, semi_global_matrix,
};
pub use types::{
    AlignError, Alignment, AlignmentMode, Column, RegionAlignError, RegionError, RegionOptions,
};

#[cfg(test)]
#[path = "../align_tests.rs"]
mod tests;
