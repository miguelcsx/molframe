//! Pairwise distance matrices and their row-block streams.

mod blocks;
mod distance;

pub use blocks::{RowBlocks, distance_rows_f32_with_context, distance_rows_with_context};
pub use distance::*;

#[cfg(test)]
mod blocks_tests;
