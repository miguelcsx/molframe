//! Pull-based batches with explicit demand and backpressure.

use super::memory::BatchGuard;
use super::{ExecutionContext, MemoryBudgetError, MemoryReservation};
use std::sync::{Arc, Mutex};

/// A bounded unit of columnar or row-oriented work.
pub trait Batch {
    /// Number of logical rows in this batch.
    fn rows(&self) -> usize;

    /// Bytes retained while this batch remains alive.
    fn retained_bytes(&self) -> usize;

    /// Reports whether this batch contains no logical rows.
    fn is_empty(&self) -> bool {
        self.rows() == 0
    }
}

/// A batch whose retained bytes are charged to an execution context.
///
/// The charge is released only when the lease is dropped, keeping accounting
/// aligned with the batch's actual lifetime.
pub struct BatchLease<B: Batch> {
    batch: LeaseBatch<B>,
    batch_guard: Option<BatchGuard>,
    reservation: Option<MemoryReservation>,
    recycler: Option<Arc<BatchBufferPool<B>>>,
}

enum LeaseBatch<B> {
    Owned(B),
    Recyclable(Arc<B>),
}

impl<B: Batch> BatchLease<B> {
    /// Charges the batch's retained bytes and creates a lease.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Exhausted`] when existing live batches and
    /// this batch would exceed the context budget.
    pub fn try_new(batch: B, context: &ExecutionContext) -> Result<Self, MemoryBudgetError> {
        let reservation = context.try_reserve(batch.retained_bytes())?;
        let batch_guard = reservation.begin_batch();
        Ok(Self {
            batch: LeaseBatch::Owned(batch),
            batch_guard: Some(batch_guard),
            reservation: Some(reservation),
            recycler: None,
        })
    }

    /// Creates a lease from capacity reserved before consuming input.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Exhausted`] if the batch retained more than
    /// the producer reserved before advancing its source.
    pub fn try_from_reservation(
        batch: B,
        mut reservation: MemoryReservation,
    ) -> Result<Self, MemoryBudgetError> {
        let retained = batch.retained_bytes();
        if retained > reservation.bytes() {
            drop(batch);
            return Err(MemoryBudgetError::Exhausted {
                requested: retained,
                available: reservation.bytes(),
            });
        }
        reservation.shrink_to(retained);
        let batch_guard = reservation.begin_batch();
        Ok(Self {
            batch: LeaseBatch::Owned(batch),
            batch_guard: Some(batch_guard),
            reservation: Some(reservation),
            recycler: None,
        })
    }

    /// Creates a lease that returns its batch capacity and memory charge to a
    /// bounded pool when dropped. The producer transfers exclusive ownership;
    /// consumers share `Arc<BatchLease<B>>` so the charge cannot be detached.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Exhausted`] if the batch retained more than
    /// the producer reserved before advancing its source.
    pub fn try_from_reservation_recycling(
        batch: B,
        mut reservation: MemoryReservation,
        recycler: Arc<BatchBufferPool<B>>,
    ) -> Result<Self, MemoryBudgetError> {
        let retained = batch.retained_bytes();
        if retained > reservation.bytes() {
            drop(batch);
            return Err(MemoryBudgetError::Exhausted {
                requested: retained,
                available: reservation.bytes(),
            });
        }
        reservation.shrink_to(retained);
        let batch_guard = reservation.begin_batch();
        Ok(Self {
            batch: LeaseBatch::Recyclable(Arc::new(batch)),
            batch_guard: Some(batch_guard),
            reservation: Some(reservation),
            recycler: Some(recycler),
        })
    }

    /// Borrows the retained batch.
    #[must_use]
    pub fn batch(&self) -> &B {
        match &self.batch {
            LeaseBatch::Owned(batch) => batch,
            LeaseBatch::Recyclable(batch) => batch,
        }
    }
}

impl<B: Batch> Drop for BatchLease<B> {
    fn drop(&mut self) {
        drop(self.batch_guard.take());
        let Some(recycler) = &self.recycler else {
            return;
        };
        let Some(reservation) = self.reservation.take() else {
            return;
        };
        let LeaseBatch::Recyclable(batch) = &self.batch else {
            return;
        };
        self.reservation = recycler.put_or_return(Arc::clone(batch), reservation);
    }
}

impl<B: Batch + std::fmt::Debug> std::fmt::Debug for BatchLease<B> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("BatchLease")
            .field("batch", self.batch())
            .field("recycling", &self.recycler.is_some())
            .finish_non_exhaustive()
    }
}

/// One-slot pool that keeps a reusable batch and its memory charge together.
///
/// The slot is intentionally singular: concurrent consumers may hold several
/// leases, but an idle source never retains capacity proportional to their
/// count. Dropping the pool releases both the allocation and its accounting.
#[derive(Debug)]
pub struct BatchBufferPool<B: Batch> {
    slot: Mutex<Option<(Arc<B>, MemoryReservation)>>,
}

impl<B: Batch> BatchBufferPool<B> {
    /// Creates an empty recycling slot.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            slot: Mutex::new(None),
        }
    }

    /// Takes exclusive ownership of the reusable batch and its live charge.
    ///
    /// Returns `None` while the previous lease is finishing destruction.
    /// No raw shared backing owner escapes the pool; consumers share the lease.
    pub fn take(&self) -> Option<(B, MemoryReservation)> {
        let Ok(mut slot) = self.slot.lock() else {
            return None;
        };
        if slot
            .as_ref()
            .is_some_and(|(batch, _)| Arc::strong_count(batch) != 1)
        {
            return None;
        }
        let (batch, reservation) = slot.take()?;
        match Arc::try_unwrap(batch) {
            Ok(batch) => Some((batch, reservation)),
            Err(batch) => {
                *slot = Some((batch, reservation));
                None
            }
        }
    }

    /// Returns an exclusively owned reusable batch and its charge. An occupied
    /// slot drops the new buffer before releasing its reservation.
    ///
    /// # Errors
    ///
    /// Rejects capacity larger than its reservation before retaining the batch.
    pub fn put(&self, batch: B, reservation: MemoryReservation) -> Result<(), MemoryBudgetError> {
        let retained = batch.retained_bytes();
        if retained > reservation.bytes() {
            drop(batch);
            return Err(MemoryBudgetError::Exhausted {
                requested: retained,
                available: reservation.bytes(),
            });
        }
        drop(self.put_or_return(Arc::new(batch), reservation));
        Ok(())
    }

    // A rejected lease still owns another Arc. Return its charge so that the
    // lease drops the backing buffer before releasing the reservation.
    fn put_or_return(
        &self,
        batch: Arc<B>,
        reservation: MemoryReservation,
    ) -> Option<MemoryReservation> {
        if let Ok(mut slot) = self.slot.lock()
            && slot.is_none()
        {
            *slot = Some((batch, reservation));
            return None;
        }
        drop(batch);
        Some(reservation)
    }
}

impl<B: Batch> Default for BatchBufferPool<B> {
    fn default() -> Self {
        Self::new()
    }
}

/// Maximum work a consumer is prepared to accept from one pull.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BatchDemand {
    /// Maximum logical rows accepted in the next batch.
    pub max_rows: usize,
    /// Maximum retained bytes accepted in the next batch.
    pub max_bytes: usize,
}

impl BatchDemand {
    /// Creates an explicit demand. Zero in either dimension means no work can
    /// currently be accepted and producers should report backpressure.
    #[must_use]
    pub const fn new(max_rows: usize, max_bytes: usize) -> Self {
        Self {
            max_rows,
            max_bytes,
        }
    }

    /// Reports whether the consumer currently has capacity for work.
    #[must_use]
    pub const fn can_accept_work(self) -> bool {
        self.max_rows > 0 && self.max_bytes > 0
    }

    /// Reports whether a batch satisfies this demand.
    #[must_use]
    pub fn accepts<B: Batch>(self, batch: &B) -> bool {
        batch.rows() <= self.max_rows && batch.retained_bytes() <= self.max_bytes
    }
}

/// Outcome of one non-blocking batch pull.
#[derive(Debug)]
pub enum Backpressure<B> {
    /// One bounded batch is ready for immediate consumption.
    Ready(B),
    /// The producer cannot currently satisfy the supplied demand.
    Pending,
    /// The producer is exhausted and will never yield another batch.
    Finished,
}

/// A pull source whose live workspace is governed by an [`ExecutionContext`].
///
/// Implementations must not return a batch exceeding `demand`. If memory is
/// temporarily unavailable, they return [`Backpressure::Pending`] without
/// blocking or allocating an unbounded queue.
pub trait BatchSource {
    /// Batch representation produced by this source.
    type Batch: Batch;
    /// Format, I/O, or computation error produced by this source.
    type Error;

    /// Attempts to produce one bounded batch.
    ///
    /// # Errors
    ///
    /// Returns the source-specific failure without consuming future input that
    /// could not be represented in the result.
    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error>;
}

#[cfg(test)]
#[path = "batch_tests.rs"]
mod tests;
