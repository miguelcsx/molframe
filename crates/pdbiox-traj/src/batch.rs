//! Budgeted pull batches over every reusable-buffer trajectory reader.

use crate::{FrameValue, Timestep, TrajectoryError, TrajectoryReader};
use pdbiox_core::execution::BatchBufferPool;
use pdbiox_core::{
    Backpressure, Batch, BatchDemand, BatchLease, BatchSource, ChunkId, DatasetId,
    ExecutionContext, LocalRow, LogicalRow, MemoryReservation,
};
use std::collections::BTreeMap;
use std::mem::size_of;
use std::sync::Arc;

// Account for both Arc control blocks and a shareable lease, even when a
// native consumer keeps its lease on the stack instead of exporting it.
const LEASE_OVERHEAD: usize = size_of::<TrajectoryBatch>() - size_of::<Timestep>()
    + size_of::<BatchLease<TrajectoryBatch>>()
    + 4 * size_of::<usize>();

/// One exactly-once frame carrying global and chunk-local identity.
#[derive(Debug)]
pub struct TrajectoryBatch {
    dataset: DatasetId,
    chunk: ChunkId,
    logical_row: LogicalRow,
    local_row: LocalRow,
    timestep: Option<Timestep>,
    retained_bytes: usize,
}

impl TrajectoryBatch {
    /// Dataset owning this frame.
    #[must_use]
    pub const fn dataset(&self) -> DatasetId {
        self.dataset
    }

    /// Stable chunk identity for this emitted frame.
    #[must_use]
    pub const fn chunk(&self) -> ChunkId {
        self.chunk
    }

    /// Stable logical row of this emitted frame.
    #[must_use]
    pub const fn logical_row(&self) -> LogicalRow {
        self.logical_row
    }

    /// Compact row inside this one-frame chunk.
    #[must_use]
    pub const fn local_row(&self) -> LocalRow {
        self.local_row
    }

    /// Borrows the reusable frame storage.
    #[must_use]
    pub fn timestep(&self) -> Option<&Timestep> {
        self.timestep.as_ref()
    }
}

impl Batch for TrajectoryBatch {
    fn rows(&self) -> usize {
        1
    }

    fn retained_bytes(&self) -> usize {
        self.retained_bytes
    }
}

// Field order keeps the reservation alive through buffer destruction on every
// early decoder, cancellation and dimension error.
struct ReservedFrame {
    timestep: Timestep,
    reservation: MemoryReservation,
}

/// Bounded `BatchSource` adapter shared by all `TrajectoryReader` formats.
#[derive(Debug)]
pub struct TrajectoryBatchSource<R> {
    reader: R,
    dataset: DatasetId,
    next_chunk: Option<ChunkId>,
    next_row: Option<LogicalRow>,
    recycler: Arc<BatchBufferPool<TrajectoryBatch>>,
    finished: bool,
}

impl<R: TrajectoryReader> TrajectoryBatchSource<R> {
    /// Wraps a reader without reading or allocating a frame.
    #[must_use]
    pub fn new(reader: R, dataset: DatasetId, first_chunk: ChunkId, first_row: LogicalRow) -> Self {
        Self {
            reader,
            dataset,
            next_chunk: Some(first_chunk),
            next_row: Some(first_row),
            recycler: Arc::new(BatchBufferPool::new()),
            finished: false,
        }
    }
}

impl<R: TrajectoryReader> BatchSource for TrajectoryBatchSource<R> {
    type Batch = TrajectoryBatch;
    type Error = TrajectoryError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error> {
        if self.finished {
            return Ok(Backpressure::Finished);
        }
        if context.cancellation().is_cancelled() {
            return Err(TrajectoryError::Cancelled);
        }
        if !demand.can_accept_work() {
            return Ok(Backpressure::Pending);
        }
        let Some(chunk) = self.next_chunk else {
            return Err(TrajectoryError::IdentityOverflow);
        };
        let Some(logical_row) = self.next_row else {
            return Err(TrajectoryError::IdentityOverflow);
        };
        let recycled = self
            .recycler
            .take()
            .filter(|(_, reservation)| reservation.belongs_to(context));
        let mut workspace = if let Some((mut batch, mut reservation)) = recycled {
            if reservation.bytes() < demand.max_bytes
                && reservation
                    .try_grow(demand.max_bytes - reservation.bytes())
                    .is_err()
            {
                self.recycler.put(batch, reservation).map_err(|_| {
                    TrajectoryError::InvalidSource {
                        format: "undercharged frame recycler",
                    }
                })?;
                return Ok(Backpressure::Pending);
            }
            let timestep = batch
                .timestep
                .take()
                .ok_or(TrajectoryError::InvalidSource {
                    format: "empty frame recycler",
                })?;
            ReservedFrame {
                timestep,
                reservation,
            }
        } else {
            let Ok(reservation) = context.try_reserve(demand.max_bytes) else {
                return Ok(Backpressure::Pending);
            };
            ReservedFrame {
                timestep: Timestep::default(),
                reservation,
            }
        };
        let decoder_bytes =
            demand
                .max_bytes
                .checked_sub(LEASE_OVERHEAD)
                .ok_or(TrajectoryError::MemoryLimit {
                    required: LEASE_OVERHEAD,
                    limit: demand.max_bytes,
                })?;
        if !self
            .reader
            .read_next_bounded(&mut workspace.timestep, decoder_bytes)?
        {
            self.finished = true;
            return Ok(Backpressure::Finished);
        }
        let retained_bytes = retained_timestep_bytes(&workspace.timestep);
        if retained_bytes > demand.max_bytes {
            return Err(TrajectoryError::MemoryLimit {
                required: retained_bytes,
                limit: demand.max_bytes,
            });
        }
        let batch = TrajectoryBatch {
            dataset: self.dataset,
            chunk,
            logical_row,
            local_row: LocalRow::new(0),
            timestep: Some(workspace.timestep),
            retained_bytes,
        };
        let lease = BatchLease::try_from_reservation_recycling(
            batch,
            workspace.reservation,
            Arc::clone(&self.recycler),
        )
        .map_err(|error| {
            let (required, limit) = match error {
                pdbiox_core::MemoryBudgetError::Exhausted {
                    requested,
                    available,
                } => (requested, available),
                pdbiox_core::MemoryBudgetError::Zero => (retained_bytes, 0),
            };
            TrajectoryError::MemoryLimit { required, limit }
        })?;
        self.next_chunk = chunk.next();
        self.next_row = logical_row.next();
        Ok(Backpressure::Ready(lease))
    }
}

fn retained_timestep_bytes(timestep: &Timestep) -> usize {
    let vectors = timestep.positions.capacity()
        + timestep.velocities.as_ref().map_or(0, Vec::capacity)
        + timestep.forces.as_ref().map_or(0, Vec::capacity);
    size_of::<Timestep>()
        .saturating_add(LEASE_OVERHEAD)
        .saturating_add(vectors.saturating_mul(size_of::<[f32; 3]>()))
        .saturating_add(retained_data_bytes(&timestep.data))
}

fn retained_data_bytes(data: &BTreeMap<Box<str>, FrameValue>) -> usize {
    data.iter().fold(0usize, |total, (key, value)| {
        let value_bytes = match value {
            FrameValue::Text(text) => text.len(),
            FrameValue::Floats(values) => values.capacity().saturating_mul(size_of::<f64>()),
            FrameValue::Float(_) | FrameValue::Integer(_) => 0,
        };
        total
            .saturating_add(size_of::<(Box<str>, FrameValue)>())
            .saturating_add(key.len())
            .saturating_add(value_bytes)
    })
}

#[cfg(test)]
#[path = "batch_tests.rs"]
mod tests;
