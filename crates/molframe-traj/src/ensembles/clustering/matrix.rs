//! Matrix validation and pair-index helpers shared by the clustering routines.

use crate::{EnsembleDistanceMatrix, EnsembleGeometryError};

/// Validates matrix dimensions, finite distances, non-negativity, and symmetry.
pub(super) fn validate_matrix(
    distances: &EnsembleDistanceMatrix,
) -> Result<(), EnsembleGeometryError> {
    let size = distances.size;
    let expected = size
        .checked_mul(size)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;

    if size == 0 || distances.values.len() != expected {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }

    for row in 0..size {
        let row_values = matrix_row(distances, row)?;
        let diagonal = row_values
            .get(row)
            .copied()
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;

        if !valid_distance(diagonal) {
            return Err(EnsembleGeometryError::InvalidParameter);
        }

        for column in row + 1..size {
            let value = row_values
                .get(column)
                .copied()
                .ok_or(EnsembleGeometryError::DimensionMismatch)?;

            let reverse = matrix_value(distances, column, row)?;

            if !valid_distance(value)
                || !valid_distance(reverse)
                || value.to_bits() != reverse.to_bits()
            {
                return Err(EnsembleGeometryError::InvalidParameter);
            }
        }
    }

    Ok(())
}

/// Returns whether one matrix entry is a valid non-negative distance.
#[inline]
fn valid_distance(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

/// Returns one complete row from the dense distance matrix.
#[inline]
pub(super) fn matrix_row(
    distances: &EnsembleDistanceMatrix,
    row: usize,
) -> Result<&[f64], EnsembleGeometryError> {
    let start = row
        .checked_mul(distances.size)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;
    let end = start
        .checked_add(distances.size)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;

    distances
        .values
        .get(start..end)
        .ok_or(EnsembleGeometryError::DimensionMismatch)
}

/// Returns one checked entry from the dense distance matrix.
#[inline]
fn matrix_value(
    distances: &EnsembleDistanceMatrix,
    row: usize,
    column: usize,
) -> Result<f64, EnsembleGeometryError> {
    let row = matrix_row(distances, row)?;

    row.get(column)
        .copied()
        .ok_or(EnsembleGeometryError::DimensionMismatch)
}

/// Returns the number of unordered pairs representable for `size` observations.
pub(super) fn checked_pair_count(size: usize) -> Option<usize> {
    size.checked_mul(size.checked_sub(1)?)?.checked_div(2)
}

/// Orders an unordered pair from the lexicographically smaller identifier.
#[inline]
pub(super) fn ordered_pair(left: usize, right: usize) -> (usize, usize) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

/// Normalizes signed zero while retaining every other floating representation.
#[inline]
pub(super) fn canonical_distance(distance: f64) -> f64 {
    if distance == 0.0 { 0.0 } else { distance }
}

/// Creates a reusable index buffer sized for one complete matrix row.
pub(super) fn reusable_index_buffer(size: usize) -> Result<Vec<usize>, EnsembleGeometryError> {
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(size)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)?;
    Ok(buffer)
}
