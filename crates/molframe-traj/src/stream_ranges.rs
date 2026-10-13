//! Parallel decode of disjoint frame ranges from seekable sources.
//!
//! Serial decode bounds sources whose frames are expensive to unpack, such as
//! compressed XTC. Here each worker owns a private reader positioned at the
//! start of one contiguous range and streams it, so decode and per-frame work
//! both scale with workers. Results reach the sink in frame order, identical
//! at every worker count, and memory is bounded by a fixed number of result
//! batches per worker rather than by the trajectory length.

use crate::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader};
use molframe_core::ExecutionContext;
use std::ops::Range;
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};

/// Results per hand-off between a worker and the consumer.
const BATCH: usize = 64;
/// Batches a worker may run ahead of the consumer.
const AHEAD: usize = 2;
/// Ranges shorter than this per worker stay on one reader.
const MIN_FRAMES_PER_WORKER: usize = 256;

type Batch<T, E> = Result<Vec<T>, E>;

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
    let share = range.len().div_ceil(workers);
    let mut total = 0_u64;
    let outcome = std::thread::scope(|scope| {
        let mut receivers = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);
        let mut first = Some(first);
        for index in 0..workers {
            let begin = range.start + index * share;
            let end = (begin + share).min(range.end);
            let (sender, receiver) = sync_channel::<Batch<T, E>>(AHEAD);
            receivers.push(receiver);
            let seeded = first.take();
            let (open, compute) = (&open, &compute);
            handles.push(scope.spawn(move || {
                let outcome = stream_range(open, seeded, begin..end, context, compute, &sender);
                if let Err(error) = outcome {
                    let _ = sender.send(Err(error));
                }
            }));
        }
        let drained = drain(&receivers, range.start, share, &mut total, &mut emit);
        drop(receivers);
        let panicked = handles.into_iter().any(|handle| handle.join().is_err());
        if panicked {
            return Err(TrajectoryError::WorkerPanicked.into());
        }
        drained
    });
    outcome?;
    if total != range.len() as u64 {
        return Err(TrajectoryError::WorkerPanicked.into());
    }
    Ok(total)
}

/// Emits batches range by range; the first error stops all later ranges.
fn drain<T, E, S>(
    receivers: &[Receiver<Batch<T, E>>],
    start: usize,
    share: usize,
    total: &mut u64,
    emit: &mut S,
) -> Result<(), E>
where
    E: From<TrajectoryError>,
    S: FnMut(usize, T) -> Result<(), E>,
{
    for (index, receiver) in receivers.iter().enumerate() {
        let mut ordinal = start + index * share;
        while let Ok(batch) = receiver.recv() {
            for value in batch? {
                emit(ordinal, value)?;
                ordinal += 1;
                *total += 1;
            }
        }
    }
    Ok(())
}

fn stream_range<R, T, E, O, C>(
    open: &O,
    seeded: Option<R>,
    range: Range<usize>,
    context: &ExecutionContext,
    compute: &C,
    sender: &SyncSender<Batch<T, E>>,
) -> Result<(), E>
where
    R: TrajectoryReader,
    E: From<TrajectoryError>,
    O: Fn() -> Result<R, TrajectoryError>,
    C: Fn(&Timestep, &ExecutionContext) -> Result<T, E>,
{
    let mut reader = if let Some(reader) = seeded {
        reader
    } else {
        let mut reader = open()?;
        reader.seek(range.start)?;
        reader
    };
    let mut frame = Timestep::default();
    let mut batch = Vec::with_capacity(BATCH);
    for _ in range.clone() {
        if context.cancellation().is_cancelled() {
            return Err(TrajectoryError::Cancelled.into());
        }
        if !reader.read_next(&mut frame)? {
            return Err(short_range().into());
        }
        batch.push(compute(&frame, context)?);
        if batch.len() == BATCH && sender.send(Ok(std::mem::take(&mut batch))).is_err() {
            return Ok(());
        }
    }
    if !batch.is_empty() {
        let _ = sender.send(Ok(batch));
    }
    Ok(())
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

const fn short_range() -> TrajectoryError {
    TrajectoryError::InvalidSource {
        format: "source ended inside the requested frame range",
    }
}

#[cfg(test)]
#[path = "stream_ranges_tests.rs"]
mod tests;
