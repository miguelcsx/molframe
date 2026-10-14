//! Blocked, deterministic, mean-centring matrix products for PCA.
//!
//! The observation matrix is one row-major `f64` buffer (`rows x cols`) that is
//! never centred in place: every product centres only the block it is working
//! on, so no second full-size copy exists and the centring itself runs on the
//! worker pool.
//!
//! Every product is cut into fixed-size blocks whose boundaries depend only on
//! the matrix shape. Blocks run on the shared worker pool, and partial sums are
//! merged in block order, so results are bit-identical at any worker count.
//!
//! Complexity, for `n` observations and `p` features:
//!
//! | product | time | extra memory per in-flight block |
//! |---|---|---|
//! | [`covariance`] | `O(n p^2)` | `p^2` partial plus one centred block |
//! | [`gram`] | `O(n^2 p)` | `n^2` partial plus one centred block |
//! | [`back_project`] | `O(n p k)` | `BLOCK x k` |
//! | [`project`] | `O(n p k)` | `BLOCK x k` |

use crate::ensemble::EnsembleGeometryError;
use faer::{Accum, Mat, MatRef, Par, Side, linalg::matmul::matmul};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{
    BlockExecutionError, BlockPlan, map_blocks_in, try_for_each_block_in,
};
use std::mem::size_of;
use std::ops::Range;

/// Observations (or features) per block.
const BLOCK: usize = 256;

/// Row-major observation matrix and the feature means subtracted on the fly.
#[derive(Clone, Copy)]
pub(crate) struct Observations<'a> {
    data: &'a [f64],
    rows: usize,
    cols: usize,
    mean: &'a [f64],
}

impl<'a> Observations<'a> {
    /// `data` must hold `rows * cols` values and `mean` must hold `cols`.
    pub(crate) fn new(
        data: &'a [f64],
        rows: usize,
        cols: usize,
        mean: &'a [f64],
    ) -> Result<Self, EnsembleGeometryError> {
        if rows.checked_mul(cols) != Some(data.len()) || mean.len() != cols {
            return Err(EnsembleGeometryError::DimensionMismatch);
        }
        Ok(Self {
            data,
            rows,
            cols,
            mean,
        })
    }

    pub(crate) const fn cols(&self) -> usize {
        self.cols
    }

    /// Centred sub-block as a freshly allocated row-major buffer.
    fn centred(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<f64> {
        let mut out = Vec::with_capacity(rows.len() * cols.len());
        for row in rows {
            let line = &self.data[row * self.cols + cols.start..row * self.cols + cols.end];
            out.extend(
                line.iter()
                    .zip(&self.mean[cols.clone()])
                    .map(|(x, m)| x - m),
            );
        }
        out
    }
}

fn block_error<E>(error: &BlockExecutionError<E>) -> EnsembleGeometryError {
    match error {
        BlockExecutionError::Memory(_) => EnsembleGeometryError::MemoryLimit,
        _ => EnsembleGeometryError::InvalidParameter,
    }
}

fn bytes(entries: usize, factor: usize) -> Result<usize, EnsembleGeometryError> {
    entries
        .checked_mul(factor)
        .and_then(|value| value.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)
}

/// `centred^T * centred / divisor`, summed over observation blocks in order.
pub(crate) fn covariance(
    observations: Observations<'_>,
    divisor: f64,
    context: &ExecutionContext,
) -> Result<Mat<f64>, EnsembleGeometryError> {
    let features = observations.cols;
    let plan = BlockPlan::new(observations.rows, BLOCK);
    let scratch = bytes(features, features.saturating_add(BLOCK))?;
    let mut total = Mat::<f64>::zeros(features, features);
    try_for_each_block_in(
        plan,
        context,
        scratch,
        |_, range| {
            let rows = range.len();
            let block = observations.centred(range, 0..features);
            let rows = MatRef::from_row_major_slice(&block, rows, features);
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
    total *= faer::Scale(divisor.recip());
    Ok(total)
}

/// `centred * centred^T / divisor`, summed over feature blocks in order.
pub(crate) fn gram(
    observations: Observations<'_>,
    divisor: f64,
    context: &ExecutionContext,
) -> Result<Mat<f64>, EnsembleGeometryError> {
    let rows_total = observations.rows;
    let plan = BlockPlan::new(observations.cols, BLOCK);
    let scratch = bytes(rows_total, rows_total.saturating_add(BLOCK))?;
    let mut total = Mat::<f64>::zeros(rows_total, rows_total);
    try_for_each_block_in(
        plan,
        context,
        scratch,
        |_, range| {
            let width = range.len();
            let block = observations.centred(0..rows_total, range);
            let columns = MatRef::from_row_major_slice(&block, rows_total, width);
            let mut partial = Mat::<f64>::zeros(rows_total, rows_total);
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
    total *= faer::Scale(divisor.recip());
    Ok(total)
}

/// `centred^T * snapshots`, one feature block per task: `p x k`.
///
/// Blocks write disjoint rows, so no reduction order is involved.
pub(crate) fn back_project(
    observations: Observations<'_>,
    snapshots: MatRef<'_, f64>,
    context: &ExecutionContext,
) -> Result<Mat<f64>, EnsembleGeometryError> {
    let count = snapshots.ncols();
    let rows_total = observations.rows;
    let plan = BlockPlan::new(observations.cols, BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range| {
        let width = range.len();
        let block = observations.centred(0..rows_total, range);
        let columns = MatRef::from_row_major_slice(&block, rows_total, width);
        let mut out = Mat::<f64>::zeros(width, count);
        matmul(
            out.as_mut(),
            Accum::Replace,
            columns.transpose(),
            snapshots,
            1.0,
            Par::Seq,
        );
        out
    })
    .map_err(|_| EnsembleGeometryError::InvalidParameter)?;
    let mut out = Mat::<f64>::zeros(observations.cols, count);
    let mut start = 0;
    for block in &blocks {
        out.as_mut()
            .subrows_mut(start, block.nrows())
            .copy_from(block);
        start += block.nrows();
    }
    Ok(out)
}

/// One row of projections per observation: `centred * components`.
///
/// `components` is `p x k`.
pub(crate) fn project(
    observations: Observations<'_>,
    components: MatRef<'_, f64>,
    context: &ExecutionContext,
) -> Result<Box<[Box<[f64]>]>, EnsembleGeometryError> {
    let count = components.ncols();
    let features = observations.cols;
    let plan = BlockPlan::new(observations.rows, BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range| {
        let rows = range.len();
        let block = observations.centred(range, 0..features);
        let source = MatRef::from_row_major_slice(&block, rows, features);
        let mut out = Mat::<f64>::zeros(rows, count);
        matmul(
            out.as_mut(),
            Accum::Replace,
            source,
            components,
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

/// Top `limit` eigenpairs of a symmetric matrix, value-descending, ties by index.
pub(crate) fn top_eigenpairs(
    matrix: &Mat<f64>,
    limit: usize,
) -> Result<(Vec<f64>, Mat<f64>), EnsembleGeometryError> {
    let decomposition = matrix
        .self_adjoint_eigen(Side::Lower)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)?;
    let values: Vec<f64> = decomposition.S().column_vector().iter().copied().collect();
    let mut indices: Vec<usize> = (0..values.len()).collect();
    indices.sort_by(|&left, &right| {
        values[right]
            .total_cmp(&values[left])
            .then(left.cmp(&right))
    });
    indices.truncate(limit);
    let vectors = Mat::from_fn(matrix.nrows(), indices.len(), |row, column| {
        decomposition.U()[(row, indices[column])]
    });
    Ok((
        indices.iter().map(|&index| values[index]).collect(),
        vectors,
    ))
}

#[cfg(test)]
#[path = "pca_linalg_tests.rs"]
mod tests;
