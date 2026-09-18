//! Incremental trajectory geometry with borrowed frames and retained outputs.

use crate::{FrameAlignment, TrajectoryError, TrajectoryReader, run_analysis_stream};
use molframe_core::{ExecutionContext, execution::Retained};

/// Streams one RMSD per frame without retaining the series.
///
/// Reference coordinates are borrowed and remain in topology order.
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
    run_analysis_stream(reader, context, frame_workspace_bytes, |frame, _| {
        let value = match alignment {
            FrameAlignment::None => molframe_geom::rmsd(&frame.positions, reference)
                .map_err(|_| TrajectoryError::DegenerateFit)?,
            FrameAlignment::Rigid => {
                molframe_geom::superpose(&frame.positions, reference)
                    .map_err(|_| TrajectoryError::DegenerateFit)?
                    .rmsd
            }
        };
        emit(frame.frame, frame.time, value)
    })
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
    let mut mean = vec![[0.0_f64; 3]; atoms];
    let mut squares = vec![0.0_f64; atoms];
    let mut seen = 0_u64;
    run_analysis_stream(reader, context, frame_workspace_bytes, |frame, _| {
        if frame.positions.len() != atoms {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: atoms,
                found: frame.positions.len(),
            });
        }
        seen += 1;
        if seen > (1_u64 << 53) {
            return Err(TrajectoryError::IdentityOverflow);
        }
        let inverse = 1.0
            / crate::numeric::f64_from_usize(
                usize::try_from(seen).map_err(|_| TrajectoryError::IdentityOverflow)?,
            )
            .ok_or(TrajectoryError::IdentityOverflow)?;
        for ((point, mean), square) in frame.positions.iter().zip(&mut mean).zip(&mut squares) {
            let delta = std::array::from_fn::<_, 3, _>(|axis| f64::from(point[axis]) - mean[axis]);
            for axis in 0..3 {
                mean[axis] += delta[axis] * inverse;
            }
            *square += (0..3)
                .map(|axis| delta[axis] * (f64::from(point[axis]) - mean[axis]))
                .sum::<f64>();
        }
        Ok(())
    })?;
    if seen == 0 {
        return Err(TrajectoryError::InvalidSource {
            format: "empty RMSF stream",
        });
    }
    let divisor = crate::numeric::f64_from_usize(
        usize::try_from(seen).map_err(|_| TrajectoryError::IdentityOverflow)?,
    )
    .ok_or(TrajectoryError::IdentityOverflow)?;
    for square in &mut squares {
        *square = (*square / divisor).max(0.0).sqrt();
    }
    drop(mean);
    let retained_bytes = squares.capacity().saturating_mul(8);
    Retained::new(squares, reservation, retained_bytes).map_err(|_| TrajectoryError::MemoryLimit {
        required: retained_bytes,
        limit: bytes,
    })
}

#[cfg(test)]
#[path = "stream_metrics_tests.rs"]
mod tests;
