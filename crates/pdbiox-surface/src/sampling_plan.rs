//! Reusable validated sampling data bound to one execution account.

use crate::{
    SasaError,
    neighbourhood::{conservative_f32, validate_radii},
    sampling::fibonacci_sphere,
};
use pdbiox_core::{ExecutionContext, MemoryReservation};
use pdbiox_spatial::SpatialError;

/// Validated radii, canonical atom identifiers and deterministic sample directions.
///
/// A trajectory reuses this data across all frames. Radii are borrowed; owned
/// indices and directions keep their reservation until the sampler is dropped.
/// Frame indices and worker scratch use the same captured execution context.
pub struct SasaSampler<'a> {
    pub(super) radii: &'a [f32],
    pub(super) probe: f64,
    pub(super) cutoff: f32,
    pub(super) targets: Vec<u32>,
    pub(super) directions: Vec<[f64; 3]>,
    pub(super) context: ExecutionContext,
    _reservation: MemoryReservation,
}

impl<'a> SasaSampler<'a> {
    /// Validates and reserves reusable sampling buffers before allocating them.
    ///
    /// # Errors
    ///
    /// Returns invalid radii, probe or samples, addressability, cancellation or
    /// shared-budget errors.
    pub fn new(
        radii: &'a [f32],
        probe: f32,
        samples: u16,
        context: &ExecutionContext,
    ) -> Result<Self, SasaError> {
        if samples == 0 {
            return Err(SasaError::NoPoints);
        }
        validate_radii(radii, probe)?;
        if context.cancellation().is_cancelled() {
            return Err(SpatialError::Cancelled.into());
        }
        let atoms = u32::try_from(radii.len()).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let samples = if atoms == 0 { 0 } else { samples };
        let bytes = radii
            .len()
            .checked_mul(size_of::<u32>())
            .and_then(|bytes| bytes.checked_add(usize::from(samples) * size_of::<[f64; 3]>()))
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let reservation = context.try_reserve(bytes).map_err(SpatialError::Memory)?;
        let targets = (0..atoms).collect();
        let directions = fibonacci_sphere(samples);
        let probe = f64::from(probe);
        let widest = radii
            .iter()
            .map(|radius| f64::from(*radius) + probe)
            .fold(0.0_f64, f64::max);
        Ok(Self {
            radii,
            probe,
            cutoff: conservative_f32(2.0 * widest),
            targets,
            directions,
            context: context.clone(),
            _reservation: reservation,
        })
    }
}

impl std::fmt::Debug for SasaSampler<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("SasaSampler")
            .field("atoms", &self.radii.len())
            .field("samples", &self.directions.len())
            .field("probe", &self.probe)
            .finish_non_exhaustive()
    }
}
