//! Blocked, deterministic matrix products for PCA.
//!
//! Every product is cut into fixed-size blocks whose boundaries depend only on
//! the matrix shape. Blocks run on the shared worker pool, and partial sums are
//! merged in block order, so results are bit-identical at any worker count.

use crate::ensemble::EnsembleGeometryError;
use faer::{Accum, Mat, MatRef, Par, linalg::matmul::matmul};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{
    BlockExecutionError, BlockPlan, map_blocks_in, try_for_each_block_in,
};
use nalgebra::DMatrix;
use std::mem::size_of;

/// Observations (or features) per block.
const BLOCK: usize = 256;

fn view(matrix: &DMatrix<f64>) -> MatRef<'_, f64> {
    MatRef::from_column_major_slice(matrix.as_slice(), matrix.nrows(), matrix.ncols())
}

fn to_nalgebra(matrix: &Mat<f64>) -> DMatrix<f64> {
    DMatrix::from_fn(matrix.nrows(), matrix.ncols(), |row, column| {
        matrix[(row, column)]
    })
}

fn block_error<E>(error: &BlockExecutionError<E>) -> EnsembleGeometryError {
    match error {
        BlockExecutionError::Memory(_) => EnsembleGeometryError::MemoryLimit,
        _ => EnsembleGeometryError::InvalidParameter,
    }
}

/// `matrix^T * matrix / divisor`, summed over observation blocks in order.
pub(crate) fn covariance(
    matrix: &DMatrix<f64>,
    divisor: f64,
    context: &ExecutionContext,
) -> Result<DMatrix<f64>, EnsembleGeometryError> {
    let features = matrix.ncols();
    let source = view(matrix);
    let plan = BlockPlan::new(matrix.nrows(), BLOCK);
    let bytes = features
        .checked_mul(features)
        .and_then(|entries| entries.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let mut total = Mat::<f64>::zeros(features, features);
    try_for_each_block_in(
        plan,
        context,
        bytes,
        |_, range| {
            let rows = source.subrows(range.start, range.len());
            let mut partial = Mat::<f64>::zeros(features, features);
            matmul(
                partial.as_mut(),
                Accum::Replace,
                rows.transpose(),
                rows,
                1.0,
                Par::Seq,
            );
            Ok::<_, EnsembleGeometryError>(partial)
        },
        |partial| {
            total += &partial;
            Ok(())
        },
    )
    .map_err(|error| block_error(&error))?;
    Ok(to_nalgebra(&total) / divisor)
}

/// `matrix * matrix^T / divisor`, summed over feature blocks in order.
pub(crate) fn gram(
    matrix: &DMatrix<f64>,
    divisor: f64,
    context: &ExecutionContext,
) -> Result<DMatrix<f64>, EnsembleGeometryError> {
    let observations = matrix.nrows();
    let source = view(matrix);
    let plan = BlockPlan::new(matrix.ncols(), BLOCK);
    let bytes = observations
        .checked_mul(observations)
        .and_then(|entries| entries.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let mut total = Mat::<f64>::zeros(observations, observations);
    try_for_each_block_in(
        plan,
        context,
        bytes,
        |_, range| {
            let columns = source.subcols(range.start, range.len());
            let mut partial = Mat::<f64>::zeros(observations, observations);
            matmul(
                partial.as_mut(),
                Accum::Replace,
                columns,
                columns.transpose(),
                1.0,
                Par::Seq,
            );
            Ok::<_, EnsembleGeometryError>(partial)
        },
        |partial| {
            total += &partial;
            Ok(())
        },
    )
    .map_err(|error| block_error(&error))?;
    Ok(to_nalgebra(&total) / divisor)
}

/// One row of projections per observation: `matrix * components^T`.
pub(crate) fn project(
    matrix: &DMatrix<f64>,
    components: &[Box<[f64]>],
    context: &ExecutionContext,
) -> Result<Box<[Box<[f64]>]>, EnsembleGeometryError> {
    let source = view(matrix);
    let count = components.len();
    let weights = Mat::<f64>::from_fn(matrix.ncols(), count, |feature, component| {
        components[component][feature]
    });
    let plan = BlockPlan::new(matrix.nrows(), BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range| {
        let mut out = Mat::<f64>::zeros(range.len(), count);
        matmul(
            out.as_mut(),
            Accum::Replace,
            source.subrows(range.start, range.len()),
            weights.as_ref(),
            1.0,
            Par::Seq,
        );
        out
    })
    .map_err(|_| EnsembleGeometryError::InvalidParameter)?;
    Ok(blocks
        .iter()
        .flat_map(|block| {
            (0..block.nrows()).map(|row| (0..count).map(|column| block[(row, column)]).collect())
        })
        .collect())
}

#[cfg(test)]
#[path = "pca_linalg_tests.rs"]
mod tests;
