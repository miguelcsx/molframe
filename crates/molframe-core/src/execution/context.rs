//! Immutable policy and shared state for one execution graph.

use super::admission::Admission;
use super::memory::MemoryAccount;
use super::spill::DiskAccount;
use super::{
    CancellationToken, MemoryBudget, MemoryBudgetError, MemoryReservation, ScratchPolicy,
    SpillError, SpillFile, TempStoragePolicy,
};
use crate::parallel::SharedPool;
use std::fmt;
use std::sync::Arc;

/// Shared controls passed through every stage of an out-of-core operation.
#[derive(Clone, Debug)]
pub struct ExecutionContext {
    pub(super) memory: Arc<MemoryAccount>,
    cancellation: CancellationToken,
    scratch: ScratchPolicy,
    temp_storage: TempStoragePolicy,
    disk: Option<Arc<DiskAccount>>,
    executor: Arc<SharedPool>,
    worker_budget: usize,
    pub(crate) admission: Arc<Admission>,
}

impl ExecutionContext {
    /// Starts a declarative context builder with bounded defaults.
    #[must_use]
    pub fn builder() -> ExecutionContextBuilder {
        ExecutionContextBuilder::default()
    }

    /// Returns this context's execution-memory ceiling.
    #[must_use]
    pub fn memory_budget(&self) -> MemoryBudget {
        self.memory.budget()
    }

    /// Maximum native tasks one operation may submit to the shared pool.
    #[must_use]
    pub const fn worker_budget(&self) -> usize {
        self.worker_budget
    }

    pub(crate) fn executor(&self) -> &SharedPool {
        &self.executor
    }

    /// Returns the shared cooperative cancellation token.
    #[must_use]
    pub const fn cancellation(&self) -> &CancellationToken {
        &self.cancellation
    }

    /// Returns the reusable scratch policy.
    #[must_use]
    pub const fn scratch_policy(&self) -> ScratchPolicy {
        self.scratch
    }

    /// Returns the temporary-storage policy.
    #[must_use]
    pub const fn temp_storage_policy(&self) -> &TempStoragePolicy {
        &self.temp_storage
    }

    /// Creates one bounded, execution-owned spill file.
    ///
    /// The key becomes the deterministic file name below the configured root.
    /// A file with the same key is never overwritten.
    ///
    /// # Errors
    ///
    /// Returns [`SpillError::Disabled`] when spill is not configured, or a
    /// storage error when the key, disk budget, or destination is invalid.
    pub fn create_spill_file(&self, key: &str) -> Result<SpillFile, SpillError> {
        let Some(account) = &self.disk else {
            return Err(SpillError::Disabled);
        };
        SpillFile::create(Arc::clone(account), key)
    }

    /// Returns temporary bytes currently retained by this execution graph.
    #[must_use]
    pub fn spill_bytes(&self) -> u64 {
        match &self.disk {
            Some(account) => account.used_bytes(),
            None => 0,
        }
    }

    /// Total physical bytes written to execution-owned spill files.
    #[must_use]
    pub fn total_spilled_bytes(&self) -> u64 {
        match &self.disk {
            Some(account) => account.spilled_bytes(),
            None => 0,
        }
    }

    /// Bytes currently retained by execution-owned allocations.
    #[must_use]
    pub fn reserved_bytes(&self) -> usize {
        self.memory.reserved_bytes()
    }

    /// Highest number of execution-owned bytes retained concurrently.
    #[must_use]
    pub fn peak_reserved_bytes(&self) -> usize {
        self.memory.peak_reserved_bytes()
    }

    /// Batch leases that are currently alive.
    #[must_use]
    pub fn live_batches(&self) -> usize {
        self.memory.live_batches()
    }

    /// Highest number of batch leases alive concurrently.
    #[must_use]
    pub fn peak_live_batches(&self) -> usize {
        self.memory.peak_live_batches()
    }

    /// Charges retained bytes against this context until the result is dropped.
    ///
    /// # Errors
    ///
    /// Returns [`MemoryBudgetError::Exhausted`] when the requested bytes do not
    /// fit alongside current reservations.
    pub fn try_reserve(&self, bytes: usize) -> Result<MemoryReservation, MemoryBudgetError> {
        MemoryAccount::reserve(&self.memory, bytes)
    }
}

impl Default for ExecutionContext {
    fn default() -> Self {
        Self::builder().build_infallible()
    }
}

/// Declarative construction of an [`ExecutionContext`].
#[derive(Clone, Debug, Default)]
pub struct ExecutionContextBuilder {
    memory_budget: MemoryBudget,
    cancellation: CancellationToken,
    scratch: ScratchPolicy,
    temp_storage: TempStoragePolicy,
    worker_budget: Option<usize>,
}

impl ExecutionContextBuilder {
    /// Replaces the default 100 MB execution budget.
    #[must_use]
    pub const fn memory_budget(mut self, budget: MemoryBudget) -> Self {
        self.memory_budget = budget;
        self
    }

    /// Uses an existing cancellation domain.
    #[must_use]
    pub fn cancellation(mut self, token: CancellationToken) -> Self {
        self.cancellation = token;
        self
    }

    /// Replaces the reusable scratch policy.
    #[must_use]
    pub const fn scratch_policy(mut self, policy: ScratchPolicy) -> Self {
        self.scratch = policy;
        self
    }

    /// Replaces the temporary-storage policy.
    #[must_use]
    pub fn temp_storage_policy(mut self, policy: TempStoragePolicy) -> Self {
        self.temp_storage = policy;
        self
    }

    /// Caps native tasks submitted by one operation to the shared worker pool.
    #[must_use]
    pub const fn worker_budget(mut self, workers: usize) -> Self {
        self.worker_budget = Some(workers);
        self
    }

    /// Validates policy relationships and constructs the context.
    ///
    /// # Errors
    ///
    /// Returns [`ContextError`] if reusable scratch could exceed the complete
    /// execution-memory budget or spill is configured with no disk capacity.
    pub fn build(self) -> Result<ExecutionContext, ContextError> {
        if self.worker_budget == Some(0) {
            return Err(ContextError::ZeroWorkers);
        }
        if self.scratch.max_bytes() > self.memory_budget.bytes() {
            return Err(ContextError::ScratchExceedsMemory {
                scratch: self.scratch.max_bytes(),
                memory: self.memory_budget.bytes(),
            });
        }
        if matches!(
            self.temp_storage,
            TempStoragePolicy::Directory { max_bytes: 0, .. }
        ) {
            return Err(ContextError::EmptyTempStorage);
        }
        Ok(self.finish())
    }

    fn build_infallible(self) -> ExecutionContext {
        self.finish()
    }

    fn finish(self) -> ExecutionContext {
        let disk = DiskAccount::from_policy(&self.temp_storage);
        let executor = SharedPool::global();
        let worker_budget = self.worker_budget.map_or_else(
            || executor.capacity(),
            |workers| workers.min(executor.capacity()),
        );
        ExecutionContext {
            memory: Arc::new(MemoryAccount::new(self.memory_budget)),
            cancellation: self.cancellation,
            scratch: self.scratch,
            temp_storage: self.temp_storage,
            disk,
            executor,
            worker_budget,
            admission: Arc::new(Admission::new(worker_budget.saturating_mul(2))),
        }
    }
}

/// Inconsistent execution policies supplied to a context builder.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContextError {
    /// Reusable scratch alone would exceed the execution-memory ceiling.
    ScratchExceedsMemory {
        /// Configured scratch bytes.
        scratch: usize,
        /// Configured complete memory budget.
        memory: usize,
    },
    /// Enabled temporary storage must have non-zero capacity.
    EmptyTempStorage,
    /// At least one worker is required to make progress.
    ZeroWorkers,
}

impl fmt::Display for ContextError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ScratchExceedsMemory { scratch, memory } => write!(
                formatter,
                "scratch policy retains {scratch} bytes but memory budget is {memory}"
            ),
            Self::EmptyTempStorage => {
                formatter.write_str("enabled temporary storage must have non-zero capacity")
            }
            Self::ZeroWorkers => formatter.write_str("execution worker budget must be non-zero"),
        }
    }
}

impl std::error::Error for ContextError {}

#[cfg(test)]
#[path = "context_tests.rs"]
mod tests;
