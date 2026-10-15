//! Read-ahead frame pipeline: serial decode overlapped with parallel per-frame work.
//!
//! A window of frames is decoded in source order on the caller's thread while
//! the previous window is computed on the shared worker pool. Results reach
//! the sink in frame order, so output is identical at any worker count and
//! live frame memory is bounded by two windows.

use crate::{Timestep, TrajectoryError, TrajectoryReader, run_analysis_stream};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use std::time::{Duration, Instant};

/// Frames per window for each worker, trading overlap against live memory.
const FRAMES_PER_WORKER: usize = 4;
const MAX_WINDOW: usize = 64;
/// Frames timed serially before choosing between serial and windowed work.
const PROBE_FRAMES: u32 = 4;
/// Mean per-frame compute time below which window hand-off costs more than it saves.
const SERIAL_BELOW: Duration = Duration::from_micros(40);

struct Fill {
    count: usize,
    error: Option<TrajectoryError>,
}

fn window_frames(context: &ExecutionContext, frame_bytes: usize) -> usize {
    let wanted = context
        .worker_budget()
        .saturating_mul(FRAMES_PER_WORKER)
        .min(MAX_WINDOW);
    let available = context
        .memory_budget()
        .bytes()
        .saturating_sub(context.reserved_bytes());
    let affordable = available / frame_bytes.max(1) / 2;
    wanted.min(affordable)
}

fn fill<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    buffers: &mut [Timestep],
    context: &ExecutionContext,
    frame_bytes: usize,
) -> Fill {
    let mut count = 0;
    for slot in buffers.iter_mut() {
        if context.cancellation().is_cancelled() {
            return Fill {
                count,
                error: Some(TrajectoryError::Cancelled),
            };
        }
        match reader.read_next_bounded(slot, frame_bytes) {
            Ok(true) => count += 1,
            Ok(false) => break,
            Err(error) => {
                return Fill {
                    count,
                    error: Some(error),
                };
            }
        }
    }
    Fill { count, error: None }
}

/// Applies `compute` to every frame and feeds results to `emit` in source order.
///
/// With a worker budget of one, or when two windows do not fit the remaining
/// memory budget, frames run serially through [`run_analysis_stream`] holding
/// one `frame_workspace_bytes` reservation. Otherwise up to
/// `2 * window * frame_workspace_bytes` is reserved for the call, where the
/// window is four frames per worker capped at 64. `compute` sees one frame at
/// a time and must be a pure function of it; `emit` runs on the caller's
/// thread. The first frames are timed serially; streams whose per-frame work
/// is too cheap to repay a hand-off stay serial. The choice never changes the
/// emitted values or their order.
///
/// # Errors
///
/// Returns decoder, resource, cancellation, worker or callback errors. Frames
/// ahead of the first failure are emitted before it is reported.
pub fn run_frame_pipeline<R, T, E, C, S>(
    reader: &mut R,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    compute: C,
    emit: S,
) -> Result<u64, E>
where
    R: TrajectoryReader + ?Sized,
    T: Send,
    E: From<TrajectoryError> + Send,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E> + Sync,
    S: FnMut(&Timestep, T) -> Result<(), E>,
{
    pipeline(
        reader,
        context,
        frame_workspace_bytes,
        SERIAL_BELOW,
        compute,
        emit,
    )
}

fn pipeline<R, T, E, C, S>(
    reader: &mut R,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    serial_below: Duration,
    compute: C,
    mut emit: S,
) -> Result<u64, E>
where
    R: TrajectoryReader + ?Sized,
    T: Send,
    E: From<TrajectoryError> + Send,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E> + Sync,
    S: FnMut(&Timestep, T) -> Result<(), E>,
{
    let window = window_frames(context, frame_workspace_bytes);
    if window < 2 || context.worker_budget() < 2 {
        return run_analysis_stream(reader, context, frame_workspace_bytes, |frame, ctx| {
            let value = compute(frame, ctx)?;
            emit(frame, value)
        });
    }
    let (mut total, probe) =
        probe_serially(reader, context, frame_workspace_bytes, &compute, &mut emit)?;
    let Some(spent) = probe else {
        return Ok(total);
    };
    if spent / PROBE_FRAMES < serial_below {
        let rest = run_analysis_stream(reader, context, frame_workspace_bytes, |frame, ctx| {
            let value = compute(frame, ctx)?;
            emit(frame, value)
        })?;
        return Ok(total + rest);
    }
    let bytes = window
        .saturating_mul(2)
        .saturating_mul(frame_workspace_bytes);
    let _reservation = context
        .try_reserve(bytes)
        .map_err(|_| TrajectoryError::MemoryLimit {
            required: bytes,
            limit: context.memory_budget().bytes(),
        })?;
    let mut current: Vec<Timestep> = (0..window).map(|_| Timestep::default()).collect();
    let mut ahead: Vec<Timestep> = (0..window).map(|_| Timestep::default()).collect();
    let mut head = fill(reader, &mut current, context, frame_workspace_bytes);
    loop {
        if head.count == 0 {
            return head.error.map_or(Ok(total), |error| Err(error.into()));
        }
        let more = head.error.is_none() && head.count == window;
        let frames = &current[..head.count];
        let (results, next) = std::thread::scope(|scope| {
            let worker = scope.spawn(|| {
                map_blocks_in(BlockPlan::new(frames.len(), 1), context, |_, range| {
                    compute(&frames[range.start], context)
                })
            });
            let next = if more {
                fill(reader, &mut ahead, context, frame_workspace_bytes)
            } else {
                Fill {
                    count: 0,
                    error: None,
                }
            };
            (worker.join(), next)
        });
        let results = results
            .map_err(|_| TrajectoryError::WorkerPanicked)?
            .map_err(|_| TrajectoryError::WorkerPanicked)?;
        for (frame, result) in frames.iter().zip(results) {
            emit(frame, result?)?;
            total = total
                .checked_add(1)
                .ok_or(TrajectoryError::IdentityOverflow)?;
        }
        if let Some(error) = head.error {
            return Err(error.into());
        }
        std::mem::swap(&mut current, &mut ahead);
        head = next;
    }
}

/// Runs the first frames serially, returning how many ran and their compute
/// time, or no time when the stream ended inside the probe.
fn probe_serially<R, T, E, C, S>(
    reader: &mut R,
    context: &ExecutionContext,
    frame_bytes: usize,
    compute: &C,
    emit: &mut S,
) -> Result<(u64, Option<Duration>), E>
where
    R: TrajectoryReader + ?Sized,
    E: From<TrajectoryError>,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E>,
    S: FnMut(&Timestep, T) -> Result<(), E>,
{
    let _reservation =
        context
            .try_reserve(frame_bytes)
            .map_err(|_| TrajectoryError::MemoryLimit {
                required: frame_bytes,
                limit: context.memory_budget().bytes(),
            })?;
    let mut frame = Timestep::default();
    let mut spent = Duration::ZERO;
    for done in 0..PROBE_FRAMES {
        if context.cancellation().is_cancelled() {
            return Err(TrajectoryError::Cancelled.into());
        }
        if !reader.read_next_bounded(&mut frame, frame_bytes)? {
            return Ok((u64::from(done), None));
        }
        let started = Instant::now();
        let value = compute(&frame, context)?;
        spent += started.elapsed();
        emit(&frame, value)?;
    }
    Ok((u64::from(PROBE_FRAMES), Some(spent)))
}

#[cfg(test)]
#[path = "stream_pipeline_tests.rs"]
mod tests;
