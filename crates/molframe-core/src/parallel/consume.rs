//! Ordered consumption of budgeted canonical blocks without collecting a dataset.

use super::ordered::ResultSlot;
use super::{BlockPlan, WorkerPanicked};
use crate::execution::admission::ExecutionDepth;
use crate::{ExecutionContext, MemoryBudgetError};
use std::fmt;
use std::mem::size_of;
use std::ops::Range;

/// A bounded block execution stopped before completing its output.
#[derive(Debug)]
pub enum BlockExecutionError<E> {
    /// No complete block fits in the remaining memory budget.
    Memory(MemoryBudgetError),
    /// Cooperative cancellation was requested.
    Cancelled,
    /// A compute worker panicked.
    Worker(WorkerPanicked),
    /// The computation or ordered consumer failed.
    Operation(E),
}

impl<E: fmt::Display> fmt::Display for BlockExecutionError<E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Memory(error) => error.fmt(formatter),
            Self::Cancelled => formatter.write_str("block execution cancelled"),
            Self::Worker(error) => error.fmt(formatter),
            Self::Operation(error) => error.fmt(formatter),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for BlockExecutionError<E> {}

/// Computes bounded windows and consumes results in canonical block order.
///
/// `bytes_per_block` covers the compute scratch and retained output together;
/// the executor additionally charges its result slots. The consumer must finish
/// using each value before returning, or transfer retained output to its own
/// reservation. No reservation is held for blocks outside the current window.
/// At most twice the useful worker count is admitted. Nested calls run inline.
///
/// # Errors
///
/// Returns cancellation, admission, worker, computation or consumer failures.
/// Previously consumed results are not rolled back; sinks publish transactionally.
pub fn try_for_each_block_in<T, E, F, C>(
    plan: BlockPlan,
    context: &ExecutionContext,
    bytes_per_block: usize,
    block: F,
    mut consume: C,
) -> Result<(), BlockExecutionError<E>>
where
    T: Send,
    E: Send,
    F: Fn(usize, Range<usize>) -> Result<T, E> + Sync,
    C: FnMut(T) -> Result<(), E>,
{
    let unit = bytes_per_block
        .checked_add(size_of::<ResultSlot<T, E>>())
        .and_then(|bytes| bytes.checked_add(size_of::<Result<T, E>>()))
        .ok_or(BlockExecutionError::Memory(MemoryBudgetError::Exhausted {
            requested: usize::MAX,
            available: 0,
        }))?
        .max(1);
    let (_depth, nested) = ExecutionDepth::enter();
    let workers = if nested {
        1
    } else {
        plan.useful_workers(context.worker_budget())
    };
    let mut start = 0;
    while start < plan.blocks() {
        if context.cancellation().is_cancelled() {
            return Err(BlockExecutionError::Cancelled);
        }
        let available = context
            .memory_budget()
            .bytes()
            .saturating_sub(context.reserved_bytes());
        let count = (available / unit)
            .min(workers.saturating_mul(2))
            .min(plan.blocks() - start);
        if count == 0 {
            return Err(BlockExecutionError::Memory(MemoryBudgetError::Exhausted {
                requested: unit,
                available,
            }));
        }
        let admission = if nested {
            None
        } else {
            Some(
                context
                    .admission
                    .acquire(count, context.cancellation())
                    .ok_or(BlockExecutionError::Cancelled)?,
            )
        };
        let count = admission.as_ref().map_or(1, |guard| guard.count);
        let mut reservation = context
            .try_reserve(count * unit)
            .map_err(BlockExecutionError::Memory)?;
        context.executor().consume_window(
            plan,
            start..start + count,
            workers,
            &block,
            |result| {
                if context.cancellation().is_cancelled() {
                    return Err(BlockExecutionError::Cancelled);
                }
                consume(result.map_err(BlockExecutionError::Operation)?)
                    .map_err(BlockExecutionError::Operation)?;
                reservation.shrink_to(reservation.bytes().saturating_sub(bytes_per_block));
                Ok(())
            },
        )?;
        drop(reservation);
        start += count;
    }
    Ok(())
}

#[cfg(test)]
#[path = "consume_tests.rs"]
mod tests;
