//! Running block plans on one process-wide bounded worker pool.
//!
//! Canonical blocks are assigned dynamically inside bounded windows. Results
//! occupy their indexed slots, so scheduling never changes reduction order.
//! Nested operations run inline instead of multiplying the shared worker pool.

use super::plan::BlockPlan;
use crate::ExecutionContext;
use rayon::prelude::*;
use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// A worker thread ended without returning its blocks.
///
/// The panic itself has already been reported by the default hook; this is the
/// value the caller receives so that a partial result is never mistaken for a
/// complete one.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct WorkerPanicked;

impl fmt::Display for WorkerPanicked {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("a worker thread panicked")
    }
}

impl std::error::Error for WorkerPanicked {}

/// Runs `block` over every block of `plan`, returning results in block order.
///
/// `workers` is the caller's, never the host's (NFR-113). One worker runs
/// serially with no thread spawned, which is what keeps a small workload from
/// paying for a thread it does not need (NFR-110).
///
/// The returned order is fixed by the plan, so it is identical at any worker
/// count. That is the whole guarantee: a caller merging these results in the
/// order they arrive is merging them in block order.
///
/// # Errors
///
/// Returns [`WorkerPanicked`] when a worker thread ends without returning.
#[cfg(test)]
pub(crate) fn map_blocks<T, F>(
    plan: BlockPlan,
    workers: usize,
    block: F,
) -> Result<Vec<T>, WorkerPanicked>
where
    T: Send,
    F: Fn(usize, std::ops::Range<usize>) -> T + Sync,
{
    SharedPool::global().map(plan, workers, block)
}

/// Runs a block plan on the pool and worker budget owned by an execution graph.
///
/// Concurrent and nested operations share the same bounded native pool. The
/// fixed block plan and final block-index ordering preserve deterministic
/// reductions independently of scheduling.
///
/// # Errors
///
/// Returns [`WorkerPanicked`] if a native task panics.
pub fn map_blocks_in<T, F>(
    plan: BlockPlan,
    context: &ExecutionContext,
    block: F,
) -> Result<Vec<T>, WorkerPanicked>
where
    T: Send,
    F: Fn(usize, std::ops::Range<usize>) -> T + Sync,
{
    context.executor().map(plan, context.worker_budget(), block)
}

pub(crate) struct SharedPool {
    pub(super) pool: Option<rayon::ThreadPool>,
    capacity: usize,
}

impl fmt::Debug for SharedPool {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SharedPool")
            .field("capacity", &self.capacity)
            .finish_non_exhaustive()
    }
}

impl SharedPool {
    pub(crate) fn global() -> Arc<Self> {
        static SHARED: OnceLock<Arc<SharedPool>> = OnceLock::new();
        Arc::clone(SHARED.get_or_init(|| Arc::new(Self::build())))
    }

    fn build() -> Self {
        let requested = std::thread::available_parallelism().map_or(1, std::num::NonZero::get);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(requested)
            .thread_name(|index| format!("pdbiox-worker-{index}"))
            .build()
            .ok();
        let capacity = pool
            .as_ref()
            .map_or(1, rayon::ThreadPool::current_num_threads);
        Self { pool, capacity }
    }

    pub(crate) const fn capacity(&self) -> usize {
        self.capacity
    }

    pub(crate) fn map<T, F>(
        &self,
        plan: BlockPlan,
        workers: usize,
        block: F,
    ) -> Result<Vec<T>, WorkerPanicked>
    where
        T: Send,
        F: Fn(usize, std::ops::Range<usize>) -> T + Sync,
    {
        let mut output = Vec::new();
        let window = plan
            .useful_workers(workers.min(self.capacity))
            .saturating_mul(2);
        for start in (0..plan.blocks()).step_by(window) {
            output.extend(self.map_window(
                plan,
                start..start.saturating_add(window).min(plan.blocks()),
                workers,
                &block,
            )?);
        }
        Ok(output)
    }

    pub(crate) fn map_window<T, F>(
        &self,
        plan: BlockPlan,
        indices: std::ops::Range<usize>,
        workers: usize,
        block: &F,
    ) -> Result<Vec<T>, WorkerPanicked>
    where
        T: Send,
        F: Fn(usize, std::ops::Range<usize>) -> T + Sync,
    {
        let workers = workers.max(1).min(self.capacity).min(indices.len());
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if workers <= 1 || self.pool.is_none() || rayon::current_thread_index().is_some() {
                return indices
                    .filter_map(|index| plan.range(index).map(|range| block(index, range)))
                    .collect();
            }
            let slots: Vec<_> = indices.clone().map(|_| Mutex::new(None)).collect();
            let next = AtomicUsize::new(indices.start);
            if let Some(pool) = &self.pool {
                pool.install(|| {
                    (0..workers).into_par_iter().for_each(|_| {
                        loop {
                            let index = next.fetch_add(1, Ordering::Relaxed);
                            if index >= indices.end {
                                break;
                            }
                            let Some(range) = plan.range(index) else {
                                break;
                            };
                            let value = block(index, range);
                            let mut slot = match slots[index - indices.start].lock() {
                                Ok(slot) => slot,
                                Err(poisoned) => poisoned.into_inner(),
                            };
                            *slot = Some(value);
                        }
                    });
                });
            }
            slots
                .into_iter()
                .filter_map(|slot| match slot.into_inner() {
                    Ok(value) => value,
                    Err(poisoned) => poisoned.into_inner(),
                })
                .collect()
        }))
        .map_err(|_| WorkerPanicked)
    }
}

#[cfg(test)]
#[path = "execute_tests.rs"]
mod tests;
