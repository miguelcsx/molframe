//! Bounded execution contracts for streaming and out-of-core work.
//!
//! Producers are pulled one batch at a time. Each pull carries explicit row and
//! byte demand, while retained memory is charged to a shared budget through an
//! RAII reservation. This keeps live workspace proportional to concurrent batch
//! demand rather than total input size. Cancellation and storage policy travel
//! in the same cheap-to-clone context.

pub(crate) mod admission;
mod batch;
mod cancellation;
mod context;
mod memory;
mod retained;
pub(crate) mod spill;
mod storage;
pub use retained::Retained;

pub use batch::BatchBufferPool;
pub use batch::{Backpressure, Batch, BatchDemand, BatchLease, BatchSource};
pub use cancellation::{CancellationToken, Cancelled};
pub use context::{ContextError, ExecutionContext, ExecutionContextBuilder};
pub use memory::{DEFAULT_MEMORY_BUDGET_BYTES, MemoryBudget, MemoryBudgetError, MemoryReservation};
pub use spill::{SpillArtifact, SpillError, SpillFile, SpillReader};
pub use storage::{ScratchPolicy, TempStoragePolicy};
