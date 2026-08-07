//! Explicit dense distance matrices.
//!
//! Constructing every pair is intrinsically quadratic, so the API names and
//! materialises that cost. Values occupy one flat row-major allocation rather
//! than one allocation per row.

use crate::distance;

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
}

/// Computes all distances within one position set.
///
/// The diagonal is exactly zero and the upper triangle is mirrored into the
/// lower triangle, avoiding duplicate calculations. Runs in `O(n²)` time and
/// space.
#[must_use]
pub fn distance_matrix(positions: &[[f32; 3]]) -> DistanceMatrix {
    let size = positions.len();
    let mut values = vec![0.0; size.saturating_mul(size)];
    for row in 0..size {
        for column in row + 1..size {
            let value = distance(positions[row], positions[column]);
            values[row * size + column] = value;
            values[column * size + row] = value;
        }
    }
    DistanceMatrix {
        rows: size,
        columns: size,
        values,
    }
}

/// Computes every distance between two position sets.
///
/// Runs in `O(left × right)` time and space. An empty input produces a matrix
/// with the corresponding zero dimension and no allocation.
#[must_use]
pub fn distance_matrix_between(left: &[[f32; 3]], right: &[[f32; 3]]) -> DistanceMatrix {
    let columns = right.len();
    let mut values = Vec::with_capacity(left.len().saturating_mul(columns));
    for left_position in left {
        for right_position in right {
            values.push(distance(*left_position, *right_position));
        }
    }
    DistanceMatrix {
        rows: left.len(),
        columns,
        values,
    }
}

#[cfg(test)]
#[path = "matrix_tests.rs"]
mod tests;
