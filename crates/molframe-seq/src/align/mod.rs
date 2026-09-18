//! Pairwise sequence alignment with affine gap costs.

mod api;
mod dynamic;
mod linear;
mod small;
mod types;

pub use api::*;
pub use types::*;

#[cfg(test)]
#[path = "../align_tests.rs"]
mod tests;
