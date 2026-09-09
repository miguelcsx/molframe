//! Caller-owned batch geometry without per-result allocation.

use std::fmt;

/// A batch geometry input or output has a different row count.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchGeometryError {
    expected: usize,
    actual: usize,
}

impl BatchGeometryError {
    /// Required row count.
    #[must_use]
    pub const fn expected(self) -> usize {
        self.expected
    }

    /// Supplied row count.
    #[must_use]
    pub const fn actual(self) -> usize {
        self.actual
    }
}

impl fmt::Display for BatchGeometryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "batch geometry expected {} rows but received {}",
            self.expected, self.actual
        )
    }
}

impl std::error::Error for BatchGeometryError {}

/// Writes pairwise distances into caller-owned storage.
///
/// Runs in `O(n)` time and `O(1)` auxiliary memory.
///
/// # Errors
///
/// Returns [`BatchGeometryError`] when an input or the output has a different
/// row count.
pub fn distances_into(
    left: &[[f32; 3]],
    right: &[[f32; 3]],
    output: &mut [f64],
) -> Result<usize, BatchGeometryError> {
    validate_lengths(left.len(), &[right.len(), output.len()])?;
    crate::simd::distances_into(left, right, output);
    Ok(left.len())
}

/// Writes angles in radians into caller-owned storage.
///
/// An undefined angle is represented by `NaN`, avoiding a per-row option or
/// validity allocation. Runs in `O(n)` time and `O(1)` auxiliary memory.
///
/// # Errors
///
/// Returns [`BatchGeometryError`] when an input or the output has a different
/// row count.
pub fn angles_into(
    first: &[[f32; 3]],
    vertices: &[[f32; 3]],
    third: &[[f32; 3]],
    output: &mut [f64],
) -> Result<usize, BatchGeometryError> {
    validate_lengths(first.len(), &[vertices.len(), third.len(), output.len()])?;
    crate::simd_measure::angles_into(first, vertices, third, output);
    Ok(first.len())
}

/// Writes signed torsions in radians into caller-owned storage.
///
/// An undefined torsion is represented by `NaN`. Runs in `O(n)` time and
/// `O(1)` auxiliary memory.
///
/// # Errors
///
/// Returns [`BatchGeometryError`] when an input or the output has a different
/// row count.
pub fn torsions_into(
    first: &[[f32; 3]],
    second: &[[f32; 3]],
    third: &[[f32; 3]],
    fourth: &[[f32; 3]],
    output: &mut [f64],
) -> Result<usize, BatchGeometryError> {
    validate_lengths(
        first.len(),
        &[second.len(), third.len(), fourth.len(), output.len()],
    )?;
    crate::simd_measure::torsions_into(first, second, third, fourth, output);
    Ok(first.len())
}

fn validate_lengths(expected: usize, lengths: &[usize]) -> Result<(), BatchGeometryError> {
    if let Some(actual) = lengths.iter().copied().find(|length| *length != expected) {
        return Err(BatchGeometryError { expected, actual });
    }
    Ok(())
}

#[cfg(test)]
#[path = "batch_measure_tests.rs"]
mod tests;
