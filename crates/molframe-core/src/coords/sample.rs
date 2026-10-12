//! Borrowed Cartesian samples evaluate one atom in O(1), without allocating a coordinate column.

use num_traits::ToPrimitive;

/// A topology-aligned linear sample of two retained coordinate columns.
/// Periodic reconstruction and alignment belong to the caller before sampling.
#[derive(Clone, Copy, Debug)]
pub struct CoordinateSample<'a> {
    start: &'a [[f32; 3]],
    end: &'a [[f32; 3]],
    weight: f32,
}

/// Invalid coordinate-sample metadata.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum CoordinateSampleError {
    /// Both columns must address the same original atom rows.
    #[error("coordinate sample requires equal-length columns")]
    Shape,
    /// Extrapolation and nonfinite interpolation weights are refused.
    #[error("coordinate sample weight must be finite and between zero and one")]
    Weight,
}

impl<'a> CoordinateSample<'a> {
    /// Validates metadata in O(1) while borrowing both original columns.
    /// Individual coordinate triples are checked when read with [`Self::position`].
    ///
    /// # Errors
    /// Refuses mismatched columns or an invalid interpolation weight.
    pub fn new(
        start: &'a [[f32; 3]],
        end: &'a [[f32; 3]],
        weight: f32,
    ) -> Result<Self, CoordinateSampleError> {
        if start.len() != end.len() {
            return Err(CoordinateSampleError::Shape);
        }
        if !weight.is_finite() || !(0.0..=1.0).contains(&weight) {
            return Err(CoordinateSampleError::Weight);
        }
        Ok(Self { start, end, weight })
    }

    /// Number of original atom addresses in this sample.
    #[must_use]
    pub const fn len(self) -> usize {
        self.start.len()
    }

    /// Whether the sample contains any original atom addresses.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.start.is_empty()
    }

    /// Reads one interpolated atom, refusing an absent row or nonfinite input.
    /// Wider intermediate arithmetic keeps opposite extreme finite positions
    /// from overflowing while preserving exact coordinates at both endpoints.
    #[must_use]
    pub fn position(self, row: usize) -> Option<[f32; 3]> {
        let start = *self.start.get(row)?;
        let end = *self.end.get(row)?;
        if !start.iter().chain(&end).all(|value| value.is_finite()) {
            return None;
        }
        if self.weight.to_bits() == 0.0f32.to_bits() || self.weight.to_bits() == (-0.0f32).to_bits()
        {
            return Some(start);
        }
        if self.weight.to_bits() == 1.0f32.to_bits() {
            return Some(end);
        }
        let weight = f64::from(self.weight);
        let mut output = [0.0; 3];
        for axis in 0..3 {
            output[axis] = (f64::from(start[axis]) * (1.0 - weight)
                + f64::from(end[axis]) * weight)
                .to_f32()?;
        }
        Some(output)
    }
}

#[cfg(test)]
#[path = "sample_tests.rs"]
mod tests;
