//! Frame-ordered total SASA over the native borrowed trajectory consumer.

use molframe_core::ExecutionContext;
use molframe_spatial::PeriodicBox;
use molframe_surface::{SasaError, SasaSampler};
use molframe_traj::{TrajectoryError, TrajectoryReader, run_analysis_stream};

/// Decoder, sampling or sink failure from a streamed surface analysis.
#[derive(Debug, thiserror::Error)]
pub enum SasaStreamError {
    /// Reading, resource admission, cancellation or sink I/O failed.
    #[error(transparent)]
    Trajectory(#[from] TrajectoryError),
    /// Surface geometry, its working set or a worker failed.
    #[error(transparent)]
    Surface(#[from] SasaError),
    /// The radius table conflicts with the policy or could not be supplied.
    #[error(transparent)]
    Radii(#[from] molframe_chem::TableError),
    /// Only Shrake-Rupley sampling can be streamed.
    #[error("streamed SASA samples points on each sphere; the lee-richards method cannot stream")]
    Method,
    /// No included atom with a known radius was left to measure.
    #[error("the topology has no included atom with a van der Waals radius to measure")]
    NoAtoms,
    /// The topology and the trajectory describe different numbers of atoms.
    #[error("the topology has {topology} atoms but the trajectory has {trajectory}")]
    TopologyAtoms {
        /// Atoms of the topology structure.
        topology: usize,
        /// Atoms of the trajectory.
        trajectory: usize,
    },
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
