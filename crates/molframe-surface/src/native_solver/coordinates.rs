//! Cartesian pose sampling while writing selected spheres, without a coordinate copy.

use super::{NativeSurfaceError, SurfaceAtoms};

/// A linear Cartesian sample between the supplied atom column and a second frame.
/// Both columns retain their original topology-aligned row addresses. Periodic
/// reconstruction or alignment must be performed upstream before supplying them.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceCoordinateSample<'a> {
    /// Borrowed coordinates of the second frame.
    pub end_positions: &'a [[f32; 3]],
    /// Interpolation weight from zero (first frame) to one (second frame).
    pub weight: f32,
}

impl SurfaceCoordinateSample<'_> {
    /// Samples one original atom address without materializing a coordinate column.
    /// Returns `None` for invalid weights, mismatched columns or nonfinite coordinates.
    #[must_use]
    pub fn position(self, start_positions: &[[f32; 3]], row: u32) -> Option<[f32; 3]> {
        molframe_core::coords::CoordinateSample::new(
            start_positions,
            self.end_positions,
            self.weight,
        )
        .ok()?
        .position(row as usize)
    }
}

#[derive(Clone, Copy)]
pub(super) struct NativeAtoms<'a> {
    pub(super) atoms: SurfaceAtoms<'a>,
    pub(super) sample: Option<SurfaceCoordinateSample<'a>>,
}

impl NativeAtoms<'_> {
    pub(super) fn validate(self) -> Result<(), NativeSurfaceError> {
        if let Some(sample) = self.sample
            && (sample.end_positions.len() != self.atoms.positions.len()
                || !sample.weight.is_finite()
                || !(0.0..=1.0).contains(&sample.weight))
        {
            return Err(NativeSurfaceError::Invalid(
                "surface sample requires aligned frames and a weight in zero to one",
            ));
        }
        Ok(())
    }

    pub(super) fn position(self, row: u32) -> Option<[f32; 3]> {
        if let Some(sample) = self.sample {
            return sample.position(self.atoms.positions, row);
        }
        let start = *self.atoms.positions.get(row as usize)?;
        if !start.iter().all(|value| value.is_finite()) {
            return None;
        }
        Some(start)
    }
}

#[cfg(test)]
#[path = "coordinates_tests.rs"]
mod tests;
