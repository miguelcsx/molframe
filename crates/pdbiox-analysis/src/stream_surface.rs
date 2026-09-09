//! Frame-ordered total SASA over the native borrowed trajectory consumer.

use pdbiox_core::ExecutionContext;
use pdbiox_spatial::PeriodicBox;
use pdbiox_surface::{SasaError, SasaSampler};
use pdbiox_traj::{TrajectoryError, TrajectoryReader, run_analysis_stream};

/// Decoder, sampling or sink failure from a streamed surface analysis.
#[derive(Debug, thiserror::Error)]
pub enum SasaStreamError {
    /// Reading, resource admission, cancellation or sink I/O failed.
    #[error(transparent)]
    Trajectory(#[from] TrajectoryError),
    /// Surface geometry, its working set or a worker failed.
    #[error(transparent)]
    Surface(#[from] SasaError),
}

/// Emits total solvent-accessible area per frame without retaining the series.
///
/// Radii are borrowed in topology order. Each frame's periodic cell is used
/// when present. Atom areas are consumed in canonical order with CPU f64
/// arithmetic, independent of worker count. No pair adjacency is materialized.
///
/// # Errors
///
/// Returns dimension, decoder, surface, resource, cancellation or sink errors.
pub fn sasa_stream<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    radii: &[f32],
    probe: f32,
    samples: u16,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    mut emit: impl FnMut(usize, Option<f64>, f64) -> Result<(), SasaStreamError>,
) -> Result<u64, SasaStreamError> {
    if radii.len() != reader.n_atoms() {
        return Err(SasaError::LengthMismatch {
            positions: reader.n_atoms(),
            radii: radii.len(),
        }
        .into());
    }
    let prepared = SasaSampler::new(radii, probe, samples, context)?;
    run_analysis_stream(reader, context, frame_workspace_bytes, |frame, _| {
        let periodic = frame
            .cell
            .map(PeriodicBox::from_cell)
            .transpose()
            .map_err(SasaError::from)?;
        let mut total = 0.0;
        prepared.visit(&frame.positions, periodic.as_ref(), |_, area| {
            total += area;
            Ok(())
        })?;
        emit(frame.frame, frame.time, total)
    })
}

#[cfg(test)]
#[path = "stream_surface_tests.rs"]
mod tests;
