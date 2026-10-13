//! Parallel decode of a frame range from seekable sources.
//!
//! Serial decode bounds sources whose frames are expensive to unpack, such as
//! compressed XTC. Each worker owns a private reader and repeatedly claims the
//! next block of frames, so decode and per-frame work scale with workers while
//! the consumer still receives results in frame order, identical at every
//! worker count. Memory is bounded by a fixed window of result blocks rather
//! than by the trajectory length.

use crate::ordered_blocks::{Ordered, plan, work};
use crate::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader};
use molframe_core::ExecutionContext;
use std::ops::Range;

/// Frames per scheduled block.
const BLOCK: usize = 64;
/// Ranges shorter than this per worker stay on one reader.
const MIN_FRAMES_PER_WORKER: usize = 256;

/// Applies `compute` to the frames in `range` and feeds results to `emit` in order.
///
/// `open` must return an independent reader over the same source each time it
/// is called. Sources without random access, ranges too short to repay extra
/// readers and worker budgets of one decode serially on a single reader. Each
/// worker holds one `frame_workspace_bytes` reservation; `compute` must be a
/// pure function of its frame and `emit` runs on the caller's thread with the
/// source frame ordinal. The range
/// is clipped to the frames the source holds.
///
/// # Errors
///
/// Returns open, seek, decoder, resource, cancellation, worker or callback
/// errors. Frames ahead of the first failure are emitted before it is reported.
pub fn run_frame_ranges<R, T, E, O, C, S>(
    open: O,
    range: Range<usize>,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    compute: C,
    mut emit: S,
) -> Result<u64, E>
where
    R: TrajectoryReader,
    T: Send,
    E: From<TrajectoryError> + Send,
    O: Fn() -> Result<R, TrajectoryError> + Sync,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E> + Sync,
    S: FnMut(usize, T) -> Result<(), E>,
{
    let mut first = open()?;
    if first.n_frames().is_some_and(|n| range.start >= n) {
        return Ok(0);
    }
    first.seek(range.start)?;
    let stop = first.n_frames().map_or(range.end, |n| range.end.min(n));
    if stop <= range.start {
        return Ok(0);
    }
    let range = range.start..stop;
    let workers = context
        .worker_budget()
        .min(range.len() / MIN_FRAMES_PER_WORKER)
        .max(1);
    if workers < 2 || first.random_access() == RandomAccess::None {
        return serial(first, range, context, frame_workspace_bytes, &compute, emit);
    }
    let bytes = frame_workspace_bytes.saturating_mul(workers);
    let _reservation = context
        .try_reserve(bytes)
        .map_err(|_| TrajectoryError::MemoryLimit {
            required: bytes,
            limit: context
                .memory_budget()
                .bytes()
                .saturating_sub(context.reserved_bytes()),
        })?;
    let mut readers = Vec::with_capacity(workers);
    readers.push(first);
    for _ in 1..workers {
        let mut reader = open()?;
        reader.seek(range.start)?;
        readers.push(reader);
    }
    let (blocks, window) = plan(range.len(), BLOCK, workers);
    let shared = Ordered::<Vec<T>, E>::new(blocks, window, workers);
    let mut total = 0_u64;
    let (drained, panicked) = std::thread::scope(|scope| {
        let handles: Vec<_> = readers
            .into_iter()
            .map(|reader| {
                let (shared, compute, range) = (&shared, &compute, &range);
                scope.spawn(move || {
                    let mut frame = Timestep::default();
                    work(shared, reader, range, BLOCK, |reader, frames| {
                        decode_block(reader, &mut frame, frames, context, compute)
                    });
                })
            })
            .collect();
        let drained = drain(&shared, range.start, &mut total, &mut emit);
        shared.stop();
        let panicked = handles.into_iter().any(|handle| handle.join().is_err());
        (drained, panicked)
    });
    if panicked {
        return Err(TrajectoryError::WorkerPanicked.into());
    }
    drained?;
    Ok(total)
}

/// Emits blocks in order; the first error stops the pipeline.
fn drain<T, E, S>(
    shared: &Ordered<Vec<T>, E>,
    start: usize,
    total: &mut u64,
    emit: &mut S,
) -> Result<(), E>
where
    E: From<TrajectoryError>,
    S: FnMut(usize, T) -> Result<(), E>,
{
    let mut ordinal = start;
    while let Some(block) = shared.next() {
        for value in block? {
            emit(ordinal, value)?;
            ordinal += 1;
            *total += 1;
        }
    }
    Ok(())
}

fn decode_block<R, T, E, C>(
    reader: &mut R,
    frame: &mut Timestep,
    frames: Range<usize>,
    context: &ExecutionContext,
    compute: &C,
) -> Result<Vec<T>, E>
where
    R: TrajectoryReader,
    E: From<TrajectoryError>,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E>,
{
    let mut out = Vec::with_capacity(frames.len());
    for _ in frames {
        if context.cancellation().is_cancelled() {
            return Err(TrajectoryError::Cancelled.into());
        }
        if !reader.read_next(frame)? {
            return Err(short_range().into());
        }
        out.push(compute(frame, context)?);
    }
    Ok(out)
}

fn serial<R, T, E, C, S>(
    mut reader: R,
    range: Range<usize>,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    compute: &C,
    mut emit: S,
) -> Result<u64, E>
where
    R: TrajectoryReader,
    E: From<TrajectoryError>,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E>,
    S: FnMut(usize, T) -> Result<(), E>,
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
    for ordinal in range {
        if context.cancellation().is_cancelled() {
            return Err(TrajectoryError::Cancelled.into());
        }
        if !reader.read_next_bounded(&mut frame, frame_workspace_bytes)? {
            return Err(short_range().into());
        }
        emit(ordinal, compute(&frame, context)?)?;
        count += 1;
    }
    Ok(count)
}

pub(crate) const fn short_range() -> TrajectoryError {
    TrajectoryError::InvalidSource {
        format: "source ended inside the requested frame range",
    }
}

#[cfg(test)]
#[path = "stream_ranges_tests.rs"]
mod tests;
