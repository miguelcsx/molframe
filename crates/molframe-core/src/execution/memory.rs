//! Shared accounting for bounded execution memory.

use std::fmt;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Default upper bound for retained execution workspace: 100 MB.
///
/// A default, not a ceiling. The caller who provisioned the machine is the only
/// party who knows what the working set may be, so a budget larger than this is
/// accepted without argument.
pub const DEFAULT_MEMORY_BUDGET_BYTES: usize = 100_000_000;

/// A validated upper bound for memory retained by one execution context.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MemoryBudget {
    bytes: usize,
}

impl MemoryBudget {
    /// Creates a non-zero budget.
    ///
    /// There is no compile-time maximum. A ceiling written into a library is a
    /// ceiling nobody can raise, and a caller reading a hundred-gigabyte entry
    /// needs a working set no constant here could have anticipated.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Zero`] when `bytes` is zero, because a
    /// zero-byte budget cannot make forward progress.
    pub const fn new(bytes: usize) -> Result<Self, MemoryBudgetError> {
        if bytes == 0 {
            return Err(MemoryBudgetError::Zero);
        }
        Ok(Self { bytes })
    }

    /// Returns the maximum number of retained bytes.
    #[must_use]
    pub const fn bytes(self) -> usize {
        self.bytes
    }
}

impl Default for MemoryBudget {
    fn default() -> Self {
        Self {
            bytes: DEFAULT_MEMORY_BUDGET_BYTES,
        }
    }
}

/// Failure to construct or reserve an execution-memory budget.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MemoryBudgetError {
    /// A zero-byte budget cannot make forward progress.
    Zero,
    /// The context has insufficient unreserved capacity.
    Exhausted {
        /// Bytes requested by the operation.
        requested: usize,
        /// Bytes still available at the time of the attempt.
        available: usize,
    },
}

impl fmt::Display for MemoryBudgetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zero => formatter.write_str("execution memory budget must be non-zero"),
            Self::Exhausted {
                requested,
                available,
            } => write!(
                formatter,
                "execution requested {requested} bytes with {available} bytes available"
            ),
        }
    }
}

impl std::error::Error for MemoryBudgetError {}

#[derive(Debug)]
pub(crate) struct MemoryAccount {
    budget: MemoryBudget,
    reserved: AtomicUsize,
    peak_reserved: AtomicUsize,
    live_batches: AtomicUsize,
    peak_live_batches: AtomicUsize,
}

impl MemoryAccount {
    pub(crate) fn new(budget: MemoryBudget) -> Self {
        Self {
            budget,
            reserved: AtomicUsize::new(0),
            peak_reserved: AtomicUsize::new(0),
            live_batches: AtomicUsize::new(0),
            peak_live_batches: AtomicUsize::new(0),
        }
    }

    pub(crate) const fn budget(&self) -> MemoryBudget {
        self.budget
    }

    pub(crate) fn reserve(
        account: &Arc<Self>,
        bytes: usize,
    ) -> Result<MemoryReservation, MemoryBudgetError> {
        let mut current = account.reserved.load(Ordering::Acquire);
        loop {
            let Some(next) = current.checked_add(bytes) else {
                return Err(MemoryBudgetError::Exhausted {
                    requested: bytes,
                    available: account.budget.bytes().saturating_sub(current),
                });
            };
            if next > account.budget.bytes() {
                return Err(MemoryBudgetError::Exhausted {
                    requested: bytes,
                    available: account.budget.bytes().saturating_sub(current),
                });
            }
            match account.reserved.compare_exchange_weak(
                current,
                next,
                Ordering::AcqRel,
                Ordering::Acquire,
            ) {
                Ok(_) => {
                    account.peak_reserved.fetch_max(next, Ordering::AcqRel);
                    return Ok(MemoryReservation {
                        account: Arc::clone(account),
                        bytes,
                    });
                }
                Err(observed) => current = observed,
            }
        }
    }

    pub(crate) fn reserved_bytes(&self) -> usize {
        self.reserved.load(Ordering::Acquire)
    }

    pub(crate) fn peak_reserved_bytes(&self) -> usize {
        self.peak_reserved.load(Ordering::Acquire)
    }

    pub(crate) fn live_batches(&self) -> usize {
        self.live_batches.load(Ordering::Acquire)
    }

    pub(crate) fn peak_live_batches(&self) -> usize {
        self.peak_live_batches.load(Ordering::Acquire)
    }

    fn begin_batch(account: &Arc<Self>) -> BatchGuard {
        let live = account.live_batches.fetch_add(1, Ordering::AcqRel) + 1;
        account.peak_live_batches.fetch_max(live, Ordering::AcqRel);
        BatchGuard {
            account: Arc::clone(account),
        }
    }
}

/// Exclusive accounting for bytes retained by an operation.
///
/// Dropping the value returns its capacity to the context. Reservations cannot
/// be cloned, preventing one charge from representing multiple owners.
#[derive(Debug)]
pub struct MemoryReservation {
    account: Arc<MemoryAccount>,
    bytes: usize,
}

impl MemoryReservation {
    /// Reports whether this reservation belongs to the supplied execution account.
    #[must_use]
    pub fn belongs_to(&self, context: &super::ExecutionContext) -> bool {
        Arc::ptr_eq(&self.account, &context.memory)
    }

    /// Returns the charged byte count.
    #[must_use]
    pub const fn bytes(&self) -> usize {
        self.bytes
    }

    /// Reduces this reservation to the exact retained byte count.
    pub fn shrink_to(&mut self, bytes: usize) {
        if bytes >= self.bytes {
            return;
        }
        let released = self.bytes - bytes;
        self.bytes = bytes;
        self.account.reserved.fetch_sub(released, Ordering::AcqRel);
    }

    /// Extends this reservation before its owner grows an allocation.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Exhausted`] without changing this
    /// reservation when the additional bytes do not fit.
    pub fn try_grow(&mut self, additional: usize) -> Result<(), MemoryBudgetError> {
        if additional == 0 {
            return Ok(());
        }
        let mut extra = MemoryAccount::reserve(&self.account, additional)?;
        let Some(bytes) = self.bytes.checked_add(additional) else {
            return Err(MemoryBudgetError::Exhausted {
                requested: additional,
                available: 0,
            });
        };
        self.bytes = bytes;
        extra.bytes = 0;
        Ok(())
    }

    pub(crate) fn begin_batch(&self) -> BatchGuard {
        MemoryAccount::begin_batch(&self.account)
    }
}

impl Drop for MemoryReservation {
    fn drop(&mut self) {
        self.account
            .reserved
            .fetch_sub(self.bytes, Ordering::AcqRel);
    }
}

/// Exactly-once accounting for one live batch lease.
#[derive(Debug)]
pub(crate) struct BatchGuard {
    account: Arc<MemoryAccount>,
}

impl Drop for BatchGuard {
    fn drop(&mut self) {
        self.account.live_batches.fetch_sub(1, Ordering::AcqRel);
    }
}

#[cfg(test)]
#[path = "memory_tests.rs"]
mod tests;
