//! Bounded ordered delivery while other canonical blocks are still computing.

use super::{BlockExecutionError, BlockPlan, SharedPool, WorkerPanicked};
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};

type Completed<T, E> = Result<Result<T, E>, WorkerPanicked>;
pub(super) type ResultSlot<T, E> = (Mutex<Option<Completed<T, E>>>, Condvar);

impl SharedPool {
    pub(super) fn consume_window<T, E, F, C>(
        &self,
        plan: BlockPlan,
        indices: Range<usize>,
        workers: usize,
        block: &F,
        mut consume: C,
    ) -> Result<(), BlockExecutionError<E>>
    where
        T: Send,
        E: Send,
        F: Fn(usize, Range<usize>) -> Result<T, E> + Sync,
        C: FnMut(Result<T, E>) -> Result<(), BlockExecutionError<E>>,
    {
        let run = |index| {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                plan.range(index).map(|range| block(index, range))
            }))
            .map_err(|_| BlockExecutionError::Worker(WorkerPanicked))
        };
        let pool = self
            .pool
            .as_ref()
            .filter(|_| workers > 1 && rayon::current_thread_index().is_none());
        let Some(pool) = pool else {
            for index in indices {
                if let Some(value) = run(index)? {
                    consume(value)?;
                }
            }
            return Ok(());
        };
        let slots: Vec<ResultSlot<T, E>> = indices
            .clone()
            .map(|_| (Mutex::new(None), Condvar::new()))
            .collect();
        let next = AtomicUsize::new(indices.start);
        // The caller stays outside the compute pool. Waiting for the next
        // canonical result cannot occupy a worker needed to produce it.
        pool.in_place_scope(|scope| {
            for _ in 0..workers.min(indices.len()) {
                scope.spawn(|_| {
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= indices.end {
                            break;
                        }
                        let value = match run(index) {
                            Ok(Some(value)) => Ok(value),
                            _ => Err(WorkerPanicked),
                        };
                        let (slot, ready) = &slots[index - indices.start];
                        let mut slot = match slot.lock() {
                            Ok(slot) => slot,
                            Err(poisoned) => poisoned.into_inner(),
                        };
                        *slot = Some(value);
                        ready.notify_one();
                    }
                });
            }
            for (slot, ready) in &slots {
                let mut slot = match slot.lock() {
                    Ok(slot) => slot,
                    Err(poisoned) => poisoned.into_inner(),
                };
                while slot.is_none() {
                    slot = match ready.wait(slot) {
                        Ok(slot) => slot,
                        Err(poisoned) => poisoned.into_inner(),
                    };
                }
                let value = slot.take();
                drop(slot);
                if let Some(value) = value {
                    consume(value.map_err(BlockExecutionError::Worker)?)?;
                }
            }
            Ok(())
        })
    }
}
