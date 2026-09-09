//! Borrowed-frame consumers with memory independent of trajectory length.

use crate::{Timestep, TrajectoryError, TrajectoryReader};
use pdbiox_core::ExecutionContext;
use std::mem::size_of;

/// Processes a trajectory in source order using one reusable frame buffer.
///
/// The reader validates each frame against the reserved workspace before
/// decoding it. The callback borrows the frame and may run parallel kernels
/// through the same execution context. Its retained state needs its own
/// reservation; borrowed frame coordinates must not escape the callback.
///
/// # Errors
///
/// Returns resource, cancellation, decoder or callback errors. A streaming sink
/// should publish its output only after this function succeeds.
pub fn run_analysis_stream<R, F, E>(
    reader: &mut R,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    mut consume: F,
) -> Result<u64, E>
where
    R: TrajectoryReader + ?Sized,
    F: FnMut(&Timestep, &ExecutionContext) -> Result<(), E>,
    E: From<TrajectoryError>,
{
    let _reservation =
        context
            .try_reserve(frame_workspace_bytes)
            .map_err(|_| TrajectoryError::MemoryLimit {
                required: frame_workspace_bytes,
                limit: context
                    .memory_budget()
                    .bytes()
                    .saturating_sub(context.reserved_bytes()),
            })?;
    let mut frame = Timestep::default();
    let mut count = 0_u64;
    loop {
        if context.cancellation().is_cancelled() {
            return Err(TrajectoryError::Cancelled.into());
        }
        if !reader.read_next_bounded(&mut frame, frame_workspace_bytes)? {
            return Ok(count);
        }
        consume(&frame, context)?;
        count = count
            .checked_add(1)
            .ok_or(TrajectoryError::IdentityOverflow)?;
    }
}

pub(crate) fn prepare_positions(
    frame: &mut Timestep,
    atoms: usize,
    bytes: usize,
) -> Result<(), TrajectoryError> {
    // Discard optional arrays before reusing a coordinates-only output buffer.
    frame.velocities = None;
    frame.forces = None;
    frame.data.clear();
    let required = atoms
        .max(frame.positions.capacity())
        .checked_mul(size_of::<[f32; 3]>())
        .and_then(|value| value.checked_add(size_of::<Timestep>()))
        .ok_or(TrajectoryError::MemoryLimit {
            required: usize::MAX,
            limit: bytes,
        })?;
    if required > bytes {
        return Err(TrajectoryError::MemoryLimit {
            required,
            limit: bytes,
        });
    }
    if atoms > frame.positions.capacity() {
        frame
            .positions
            .try_reserve_exact(atoms.saturating_sub(frame.positions.len()))
            .map_err(|_| TrajectoryError::MemoryLimit {
                required,
                limit: bytes,
            })?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "stream_consume_tests.rs"]
mod tests;
