//! Explicit dense distance matrices.
//!
//! Constructing every pair is intrinsically quadratic, so the API names and
//! materialises that cost. Values occupy one flat row-major allocation rather
//! than one allocation per row.

use crate::distance;
use std::fmt;

/// Why a dense distance matrix could not be materialised.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixError {
    /// The requested dimensions cannot be represented by `usize`.
    SizeOverflow,
    /// The contiguous matrix allocation was refused.
    AllocationFailed,
}

impl fmt::Display for MatrixError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeOverflow => formatter.write_str("distance matrix dimensions overflow usize"),
            Self::AllocationFailed => formatter.write_str("distance matrix allocation was refused"),
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
/// The diagonal is exactly zero and the upper triangle is mirrored into the
/// lower triangle, avoiding duplicate calculations. Runs in `O(n²)` time and
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
    for row in 0..size {
        for column in row + 1..size {
            let value = distance(positions[row], positions[column]);
            values[row * size + column] = value;
            values[column * size + row] = value;
        }
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
    for left_position in left {
        for right_position in right {
            values.push(distance(*left_position, *right_position));
        }
    }
    Ok(DistanceMatrix {
        rows: left.len(),
        columns,
        values,
    })
}

#[cfg(test)]
#[path = "distance_tests.rs"]
mod tests;
