//! Explicit materialization of bounded sampled atom areas.

use crate::{SasaError, SasaSampler};
use molframe_core::{ExecutionContext, execution::Retained};
use molframe_spatial::{PeriodicBox, SpatialError};

/// Collects areas whose output allocation stays charged until its final owner drops.
///
/// Input coordinates and radii are borrowed. The bounded sampler never stores
/// pair adjacency; use `visit_shrake_rupley` when an output array is unnecessary.
///
/// # Errors
///
/// Returns input, memory, cancellation, geometry or worker errors.
pub fn collect_shrake_rupley(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    samples: u16,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Retained<Vec<f64>>, SasaError> {
    SasaSampler::new(radii, probe, samples, context)?.collect(positions, periodic)
}

impl SasaSampler<'_> {
    /// Materializes one frame's areas while retaining the output reservation.
    ///
    /// # Errors
    ///
    /// Returns the same failures as `visit`, including output admission failures.
    pub fn collect(
        &self,
        positions: &[[f32; 3]],
        periodic: Option<&PeriodicBox>,
    ) -> Result<Retained<Vec<f64>>, SasaError> {
        if positions.len() != self.radii.len() {
            return Err(SasaError::LengthMismatch {
                positions: positions.len(),
                radii: self.radii.len(),
            });
        }
        let bytes = positions
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let reservation = self
            .context
            .try_reserve(bytes)
            .map_err(SpatialError::Memory)?;
        let mut areas = Vec::with_capacity(positions.len());
        self.visit(positions, periodic, |_, area| {
            areas.push(area);
            Ok(())
        })?;
        Retained::new(areas, reservation, bytes).map_err(|error| SpatialError::Memory(error).into())
    }
}
