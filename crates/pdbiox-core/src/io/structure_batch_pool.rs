//! Reservation-preserving reuse for one structural batch source.

use super::{StructureBatch, StructureBatchBuffer, StructureBatchBuilder, StructureBatchError};
use crate::execution::BatchBufferPool;
use crate::{Batch, BatchLease, ChunkDescriptor, ExecutionContext, MemoryReservation};
use std::sync::Arc;

/// One reusable columnar batch owned by a structural source.
///
/// Retained capacity never becomes invisible to the execution budget: the
/// pool stores the live reservation beside the buffer and releases both when
/// the source is dropped.
#[derive(Clone, Debug)]
pub struct StructureBatchPool {
    inner: Arc<BatchBufferPool<StructureBatch>>,
}

impl StructureBatchPool {
    /// Creates an empty one-slot pool.
    #[must_use]
    pub fn new() -> Self {
        Self {
            inner: Arc::new(BatchBufferPool::new()),
        }
    }

    /// Acquires a fresh or recycled builder after charging its full working
    /// allowance. A recycled allocation may exceed `reserve_bytes` when its
    /// existing reservation still covers it and it remains within
    /// `max_retained_bytes`. `None` means the context cannot currently make
    /// progress.
    ///
    /// # Errors
    ///
    /// Rejects unrepresentable batch identity.
    pub fn acquire(
        &self,
        descriptor: ChunkDescriptor,
        reserve_bytes: usize,
        max_retained_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Option<(StructureBatchBuilder, MemoryReservation)>, StructureBatchError> {
        let capacity = descriptor.rows();
        let recycled = self.inner.take().filter(|(batch, reservation)| {
            reservation.belongs_to(context) && batch.can_reuse(capacity, max_retained_bytes)
        });
        if let Some((batch, mut reservation)) = recycled {
            if reservation.bytes() < reserve_bytes
                && reservation
                    .try_grow(reserve_bytes - reservation.bytes())
                    .is_err()
            {
                self.inner.put(batch, reservation)?;
                return Ok(None);
            }
            let buffer = StructureBatchBuffer::new(batch);
            let builder = StructureBatchBuilder::reuse(
                buffer,
                descriptor.dataset(),
                descriptor.chunk(),
                descriptor.logical_start(),
            )?;
            return Ok(Some((builder, reservation)));
        }
        let Ok(reservation) = context.try_reserve(reserve_bytes) else {
            return Ok(None);
        };
        let builder = StructureBatchBuilder::new(
            descriptor.dataset(),
            descriptor.chunk(),
            descriptor.logical_start(),
            capacity,
        )?;
        Ok(Some((builder, reservation)))
    }

    /// Publishes a batch whose capacity returns to this pool with its lease.
    ///
    /// # Errors
    ///
    /// Rejects a batch larger than the producer's reservation.
    pub fn lease(
        &self,
        buffer: StructureBatchBuffer,
        reservation: MemoryReservation,
    ) -> Result<BatchLease<StructureBatch>, StructureBatchError> {
        BatchLease::try_from_reservation_recycling(
            buffer.into_owned(),
            reservation,
            Arc::clone(&self.inner),
        )
        .map_err(StructureBatchError::from)
    }

    /// Keeps a rejected candidate available for a subsequent demand.
    pub fn recycle(&self, buffer: StructureBatchBuffer, reservation: MemoryReservation) {
        if buffer.batch().retained_bytes() > reservation.bytes() {
            drop(buffer);
            return;
        }
        // Capacity was checked above. Rejection still destroys the batch before
        // releasing its charge, so no unaccounted idle allocation survives.
        let _ = self.inner.put(buffer.into_owned(), reservation);
    }
}

impl Default for StructureBatchPool {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[path = "structure_batch_pool_tests.rs"]
mod tests;
