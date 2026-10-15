//! Incremental trajectory geometry with borrowed frames and retained outputs.

use crate::{
    FrameAlignment, TrajectoryError, TrajectoryReader, run_analysis_stream, run_frame_pipeline,
};
use molframe_core::{ExecutionContext, execution::Retained};

/// Streams one RMSD per frame without retaining the series.
///
/// Reference coordinates are borrowed and remain in topology order. Frames
/// are fitted in parallel behind a read-ahead window when the context allows
/// more than one worker; values and their order do not depend on the count.
///
/// # Errors
///
/// Returns decoding, memory, dimension, fitting or sink errors.
pub fn rmsd_stream<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    reference: &[[f32; 3]],
    alignment: FrameAlignment,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    mut emit: impl FnMut(usize, Option<f64>, f64) -> Result<(), TrajectoryError>,
) -> Result<u64, TrajectoryError> {
    if reference.len() != reader.n_atoms() {
        return Err(TrajectoryError::AtomCountMismatch {
            expected: reader.n_atoms(),
            found: reference.len(),
        });
    }
    run_frame_pipeline(
        reader,
        context,
        frame_workspace_bytes,
        |frame, _| match alignment {
            FrameAlignment::None => molframe_geom::rmsd(&frame.positions, reference)
                .map_err(|_| TrajectoryError::DegenerateFit),
            FrameAlignment::Rigid => molframe_geom::rmsd_after_fit(&frame.positions, reference)
                .map_err(|_| TrajectoryError::DegenerateFit),
        },
        |frame, value| emit(frame.frame, frame.time, value),
    )
}

/// Computes per-atom RMSF in one pass using Welford's recurrence.
///
/// Working memory is `O(atoms)` independently of the number of frames. No
/// alignment is performed implicitly. The returned result keeps its allocation
/// charged until the final owner drops it.
///
/// # Errors
///
/// Returns resource, decoder, empty-stream or dimension errors.
pub fn rmsf_stream<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
) -> Result<Retained<Vec<f64>>, TrajectoryError> {
    let atoms = reader.n_atoms();
    let bytes = atoms.checked_mul(32).ok_or(TrajectoryError::MemoryLimit {
        required: usize::MAX,
        limit: context.memory_budget().bytes(),
    })?;
    let reservation = context
        .try_reserve(bytes)
        .map_err(|_| TrajectoryError::MemoryLimit {
            required: bytes,
            limit: context.memory_budget().bytes(),
        })?;
    let mut accumulator = molframe_geom::RmsfAccumulator::new(atoms);
    run_analysis_stream(reader, context, frame_workspace_bytes, |frame, _| {
        if frame.positions.len() != atoms {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: atoms,
                found: frame.positions.len(),
            });
        }
        accumulator
            .push(&frame.positions)
            .map_err(|_| TrajectoryError::IdentityOverflow)
    })?;
    let values = accumulator
        .finish()
        .map_err(|_| TrajectoryError::InvalidSource {
            format: "empty RMSF stream",
        })?;
    let retained_bytes = values.capacity().saturating_mul(8);
    Retained::new(values, reservation, retained_bytes).map_err(|_| TrajectoryError::MemoryLimit {
        required: retained_bytes,
        limit: bytes,
    })
}

#[cfg(test)]
#[path = "stream_metrics_tests.rs"]
mod tests;
