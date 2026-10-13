//! Ordered dense chunks decoded by background workers.
//!
//! An owned iterator over `DenseChunk` blocks of a frame range. Worker threads
//! each hold a private reader and decode interleaved chunks while the caller
//! consumes them in frame order, so a slow consumer stalls the workers once
//! the read-ahead window is full instead of buffering the trajectory. The
//! window's memory is reserved from the execution context up front, workers
//! stop when the iterator is dropped or the context is cancelled, and chunk
//! contents are identical at every worker count.

use crate::ordered_blocks::{Ordered, plan, work};
use crate::read_dense_chunk;
use crate::stream_ranges::{clip_range, short_range};
use crate::{DenseChunk, RandomAccess, Timestep, TrajectoryError, TrajectoryReader};
use molframe_core::{ExecutionContext, MemoryReservation};
use std::ops::Range;
use std::sync::Arc;
use std::thread::JoinHandle;

type Pipe = Ordered<DenseChunk, TrajectoryError>;

/// In-order iterator over chunks decoded in parallel.
pub struct ParallelChunks {
    shared: Arc<Pipe>,
    handles: Vec<JoinHandle<()>>,
    _reservation: Option<MemoryReservation>,
}

impl ParallelChunks {
    /// Starts workers over `range`, clipped to the frames the source holds.
    ///
    /// `open` must return an independent reader over the same source on each
    /// call. Sources without random access use one worker. Memory for the
    /// read-ahead window is reserved from `context`, which also supplies the
    /// worker budget and cancellation.
    ///
    /// # Errors
    ///
    /// Returns open or seek errors, or [`TrajectoryError::MemoryLimit`] when
    /// the window does not fit the memory budget.
    pub fn open<R, O>(
        open: O,
        range: Range<usize>,
        chunk_frames: usize,
        context: &ExecutionContext,
    ) -> Result<Self, TrajectoryError>
    where
        R: TrajectoryReader + Send + 'static,
        O: Fn() -> Result<R, TrajectoryError>,
    {
        let chunk_frames = chunk_frames.max(1);
        let mut first = open()?;
        let Some(range) = clip_range(&mut first, range)? else {
            return Ok(Self::empty());
        };
        let blocks = range.len().div_ceil(chunk_frames);
        let workers = if first.random_access() == RandomAccess::None {
            1
        } else {
            context.worker_budget().min(blocks).max(1)
        };
        let (blocks, window) = plan(range.len(), chunk_frames, workers);
        let chunk_bytes =
            chunk_frames.saturating_mul(first.n_atoms().saturating_mul(12).saturating_add(64));
        let bytes = chunk_bytes.saturating_mul(window + workers);
        let reservation = context
            .try_reserve(bytes)
            .map_err(|_| TrajectoryError::MemoryLimit {
                required: bytes,
                limit: context
                    .memory_budget()
                    .bytes()
                    .saturating_sub(context.reserved_bytes()),
            })?;
        let mut readers = vec![first];
        for _ in 1..workers {
            let mut reader = open()?;
            reader.seek(range.start)?;
            readers.push(reader);
        }
        let shared = Arc::new(Pipe::new(blocks, window, workers));
        let range = Arc::new(range);
        let handles = readers
            .into_iter()
            .map(|reader| {
                let (shared, range, context) =
                    (Arc::clone(&shared), Arc::clone(&range), context.clone());
                std::thread::spawn(move || {
                    let mut frame = Timestep::default();
                    work(&shared, reader, &range, chunk_frames, |reader, frames| {
                        decode_chunk(reader, &mut frame, frames, &context)
                    });
                })
            })
            .collect();
        Ok(Self {
            shared,
            handles,
            _reservation: Some(reservation),
        })
    }

    fn empty() -> Self {
        Self {
            shared: Arc::new(Pipe::new(0, 1, 0)),
            handles: Vec::new(),
            _reservation: None,
        }
    }
}

fn decode_chunk<R: TrajectoryReader>(
    reader: &mut R,
    frame: &mut Timestep,
    frames: Range<usize>,
    context: &ExecutionContext,
) -> Result<DenseChunk, TrajectoryError> {
    if context.cancellation().is_cancelled() {
        return Err(TrajectoryError::Cancelled);
    }
    let chunk = read_dense_chunk(reader, frame, frames.len())?;
    if chunk.frames() == frames.len() {
        Ok(chunk)
    } else {
        Err(short_range())
    }
}

impl std::fmt::Debug for ParallelChunks {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ParallelChunks")
            .field("workers", &self.handles.len())
            .finish_non_exhaustive()
    }
}

impl Iterator for ParallelChunks {
    type Item = Result<DenseChunk, TrajectoryError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.shared.next()
    }
}

impl Drop for ParallelChunks {
    fn drop(&mut self) {
        self.shared.stop();
        for handle in self.handles.drain(..) {
            let _ = handle.join();
        }
    }
}

#[cfg(test)]
#[path = "stream_chunks_tests.rs"]
mod tests;
