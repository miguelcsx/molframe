//! Explicit dense distance matrices.
//!
//! Constructing every pair is intrinsically quadratic, so the API names and
//! materialises that cost. Values occupy one flat row-major allocation rather
//! than one allocation per row.

use crate::distance;
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use std::fmt;
use wide::f64x4;

const PARALLEL_MATRIX_ROWS: usize = 256;
const ROWS_PER_BLOCK: usize = 32;

/// Why a dense distance matrix could not be materialised.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixError {
    /// The requested dimensions cannot be represented by `usize`.
    SizeOverflow,
    /// The contiguous matrix allocation was refused.
    AllocationFailed,
    /// A caller-supplied output slice has the wrong length.
    OutputLength,
    /// The execution context refused the estimated peak allocation.
    MemoryBudget,
    /// Cooperative cancellation was requested.
    Cancelled,
    /// A worker panicked while computing a row block.
    WorkerPanicked,
}

impl fmt::Display for MatrixError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeOverflow => formatter.write_str("distance matrix dimensions overflow usize"),
            Self::AllocationFailed => formatter.write_str("distance matrix allocation was refused"),
            Self::OutputLength => {
                formatter.write_str("distance matrix output has the wrong length")
            }
            Self::MemoryBudget => {
                formatter.write_str("distance matrix exceeds the execution memory budget")
            }
            Self::Cancelled => formatter.write_str("distance matrix execution was cancelled"),
            Self::WorkerPanicked => formatter.write_str("distance matrix worker panicked"),
        }
    }
}

impl std::error::Error for MatrixError {}

/// A dense row-major matrix of distances in ångström.
#[derive(Clone, Debug, PartialEq)]
pub struct DistanceMatrix {
    rows: usize,
    columns: usize,
    values: Vec<f64>,
}

impl DistanceMatrix {
    /// Number of rows.
    #[must_use]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    #[must_use]
    pub const fn columns(&self) -> usize {
        self.columns
    }

    /// One matrix entry, or `None` when either index is outside the matrix.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        let position = row.checked_mul(self.columns)?.checked_add(column)?;
        self.values.get(position).copied()
    }

    /// The contiguous row-major values.
    #[must_use]
    pub fn as_slice(&self) -> &[f64] {
        &self.values
    }

    /// Consumes the matrix and returns its contiguous row-major values.
    ///
    /// This is useful for FFI boundaries that can adopt the allocation without
    /// cloning its `O(rows × columns)` payload.
    #[must_use]
    pub fn into_values(self) -> Vec<f64> {
        self.values
    }
}

/// Computes all distances within one position set.
///
/// The diagonal is exactly zero. Rows are evaluated independently so this
/// routine and [`distance_matrix_with_context`] share the same vectorized hot
/// path and deterministic row-block decomposition. Runs in `O(n²)` time and
/// space.
///
/// # Errors
///
/// Returns [`MatrixError::SizeOverflow`] when the square element count cannot
/// be represented, or [`MatrixError::AllocationFailed`] when its single
/// contiguous allocation is refused.
pub fn distance_matrix(positions: &[[f32; 3]]) -> Result<DistanceMatrix, MatrixError> {
    let size = positions.len();
    let count = size.checked_mul(size).ok_or(MatrixError::SizeOverflow)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| MatrixError::AllocationFailed)?;
    values.resize(count, 0.0);
    distance_matrix_into(positions, &mut values)?;
    Ok(DistanceMatrix {
        rows: size,
        columns: size,
        values,
    })
}

/// Computes a square distance matrix into caller-owned row-major storage.
///
/// The caller can recycle `output` across frames and avoid all allocation.
/// Cache-sized row blocks and four-lane inner batches are used without changing
/// the deterministic row-major result.
///
/// # Errors
///
/// Returns [`MatrixError::SizeOverflow`] or [`MatrixError::OutputLength`].
pub fn distance_matrix_into(positions: &[[f32; 3]], output: &mut [f64]) -> Result<(), MatrixError> {
    let size = positions.len();
    let count = size.checked_mul(size).ok_or(MatrixError::SizeOverflow)?;
    if output.len() != count {
        return Err(MatrixError::OutputLength);
    }
    for (row, values) in output.chunks_exact_mut(size).enumerate() {
        distances_from(positions[row], positions, values);
    }
    Ok(())
}

/// Materialises a dense matrix under an explicit execution context.
///
/// Small inputs stay on the serial allocation-free-into path. Large inputs
/// use fixed canonical row blocks on the shared bounded worker pool; block
/// boundaries and assembly order do not depend on worker count.
///
/// # Errors
///
/// Returns a matrix, cancellation, worker or memory-admission error.
pub fn distance_matrix_with_context(
    positions: &[[f32; 3]],
    context: &ExecutionContext,
) -> Result<DistanceMatrix, MatrixError> {
    let size = positions.len();
    let count = size.checked_mul(size).ok_or(MatrixError::SizeOverflow)?;
    let bytes = count
        .checked_mul(std::mem::size_of::<f64>())
        .and_then(|value| value.checked_mul(2))
        .ok_or(MatrixError::SizeOverflow)?;
    let _reservation = context
        .try_reserve(bytes)
        .map_err(|_| MatrixError::MemoryBudget)?;
    if context.cancellation().is_cancelled() {
        return Err(MatrixError::Cancelled);
    }
    if size < PARALLEL_MATRIX_ROWS || context.worker_budget() == 1 {
        return distance_matrix(positions);
    }
    let blocks = map_blocks_in(BlockPlan::new(size, ROWS_PER_BLOCK), context, |_, rows| {
        let mut values = vec![0.0; rows.len() * size];
        for (local, row) in rows.enumerate() {
            distances_from(
                positions[row],
                positions,
                &mut values[local * size..(local + 1) * size],
            );
        }
        values
    })
    .map_err(|_| MatrixError::WorkerPanicked)?;
    if context.cancellation().is_cancelled() {
        return Err(MatrixError::Cancelled);
    }
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| MatrixError::AllocationFailed)?;
    for mut block in blocks {
        values.append(&mut block);
    }
    Ok(DistanceMatrix {
        rows: size,
        columns: size,
        values,
    })
}

/// Computes every distance between two position sets.
///
/// Runs in `O(left × right)` time and space. An empty input produces a matrix
/// with the corresponding zero dimension and no allocation.
///
/// # Errors
///
/// Returns [`MatrixError::SizeOverflow`] when the rectangular element count
/// cannot be represented, or [`MatrixError::AllocationFailed`] when its single
/// contiguous allocation is refused.
pub fn distance_matrix_between(
    left: &[[f32; 3]],
    right: &[[f32; 3]],
) -> Result<DistanceMatrix, MatrixError> {
    let columns = right.len();
    let count = left
        .len()
        .checked_mul(columns)
        .ok_or(MatrixError::SizeOverflow)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| MatrixError::AllocationFailed)?;
    values.resize(count, 0.0);
    for (left_position, row) in left.iter().zip(values.chunks_exact_mut(columns)) {
        distances_from(*left_position, right, row);
    }
    Ok(DistanceMatrix {
        rows: left.len(),
        columns,
        values,
    })
}

fn distances_from(origin: [f32; 3], positions: &[[f32; 3]], output: &mut [f64]) {
    let [origin_x, origin_y, origin_z] = origin.map(f64::from);
    let (chunks, remainder) = positions.as_chunks::<4>();
    let (output_chunks, output_remainder) = output.as_chunks_mut::<4>();
    for (points, values) in chunks.iter().zip(output_chunks) {
        let x = f64x4::new([
            f64::from(points[0][0]),
            f64::from(points[1][0]),
            f64::from(points[2][0]),
            f64::from(points[3][0]),
        ]);
        let y = f64x4::new([
            f64::from(points[0][1]),
            f64::from(points[1][1]),
            f64::from(points[2][1]),
            f64::from(points[3][1]),
        ]);
        let z = f64x4::new([
            f64::from(points[0][2]),
            f64::from(points[1][2]),
            f64::from(points[2][2]),
            f64::from(points[3][2]),
        ]);
        let dx = x - f64x4::splat(origin_x);
        let dy = y - f64x4::splat(origin_y);
        let dz = z - f64x4::splat(origin_z);
        values.copy_from_slice((dx * dx + dy * dy + dz * dz).sqrt().as_array());
    }
    for (point, value) in remainder.iter().zip(output_remainder) {
        *value = distance(origin, *point);
    }
}

#[cfg(test)]
#[path = "distance_tests.rs"]
mod tests;
