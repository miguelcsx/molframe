//! Pairwise distance matrices and their row-block streams.

mod blocks;
mod distance;

pub use blocks::{RowBlocks, distance_rows_f32_with_context, distance_rows_with_context};
pub use distance::{
    DistanceMatrix, MatrixError, condensed_distances_f32_with_context,
    condensed_distances_with_context, distance_block_rows, distance_matrix,
    distance_matrix_between, distance_matrix_f32_with_context, distance_matrix_into,
    distance_matrix_with_context, distance_rows_f32_into, distance_rows_into,
};

#[cfg(test)]
mod blocks_tests;
