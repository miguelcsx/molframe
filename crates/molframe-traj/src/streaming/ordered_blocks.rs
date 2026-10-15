//! Ordered hand-off of frame blocks that workers decode out of order.
//!
//! Blocks are numbered in frame order and claimed from one shared counter, so
//! every worker always holds the next block nobody has taken, and decode of
//! consecutive blocks overlaps across workers. A claim is refused while the
//! block would sit more than `window` blocks ahead of the consumer, which
//! bounds the results in flight, and the consumer receives blocks strictly in
//! order whatever the worker count, so output is deterministic.

use crate::{TrajectoryError, TrajectoryReader};
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::{Condvar, Mutex, MutexGuard};

struct State<T, E> {
    claimed: usize,
    emitted: usize,
    live: usize,
    stopped: bool,
    ready: BTreeMap<usize, Result<T, E>>,
}

/// Shared scheduler between decode workers and one in-order consumer.
pub(crate) struct Ordered<T, E> {
    state: Mutex<State<T, E>>,
    wake: Condvar,
    blocks: usize,
    window: usize,
}

impl<T, E: From<TrajectoryError>> Ordered<T, E> {
    pub(crate) fn new(blocks: usize, window: usize, workers: usize) -> Self {
        Self {
            state: Mutex::new(State {
                claimed: 0,
                emitted: 0,
                live: workers,
                stopped: false,
                ready: BTreeMap::new(),
            }),
            wake: Condvar::new(),
            blocks,
            window: window.max(1),
        }
    }

    fn lock(&self) -> MutexGuard<'_, State<T, E>> {
        match self.state.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        }
    }

    fn wait<'a>(&self, guard: MutexGuard<'a, State<T, E>>) -> MutexGuard<'a, State<T, E>> {
        match self.wake.wait(guard) {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        }
    }

    /// Takes the next unclaimed block, waiting while the window is full.
    fn claim(&self) -> Option<usize> {
        let mut state = self.lock();
        loop {
            if state.stopped || state.claimed >= self.blocks {
                return None;
            }
            if state.claimed < state.emitted + self.window {
                state.claimed += 1;
                return Some(state.claimed - 1);
            }
            state = self.wait(state);
        }
    }

    fn complete(&self, block: usize, result: Result<T, E>) {
        self.lock().ready.insert(block, result);
        self.wake.notify_all();
    }

    fn retire(&self) {
        self.lock().live -= 1;
        self.wake.notify_all();
    }

    /// Releases every worker; later claims return nothing.
    pub(crate) fn stop(&self) {
        self.lock().stopped = true;
        self.wake.notify_all();
    }

    /// The next block in frame order, or `None` once all were delivered.
    pub(crate) fn next(&self) -> Option<Result<T, E>> {
        let mut state = self.lock();
        loop {
            if state.emitted >= self.blocks {
                return None;
            }
            let wanted = state.emitted;
            if let Some(result) = state.ready.remove(&wanted) {
                state.emitted += 1;
                if result.is_err() {
                    state.emitted = self.blocks;
                }
                self.wake.notify_all();
                return Some(result);
            }
            if state.live == 0 || state.stopped {
                state.emitted = self.blocks;
                return Some(Err(TrajectoryError::WorkerPanicked.into()));
            }
            state = self.wait(state);
        }
    }
}

struct Retire<'a, T, E: From<TrajectoryError>>(&'a Ordered<T, E>);

impl<T, E: From<TrajectoryError>> Drop for Retire<'_, T, E> {
    fn drop(&mut self) {
        self.0.retire();
    }
}

/// Decodes claimed blocks of `frames` with `reader` until none remain.
///
/// `decode` turns one block's frame range into a result and may assume the
/// reader is positioned at its start. A failure is delivered in its block's
/// position and ends this worker.
pub(crate) fn work<R, T, E, F>(
    shared: &Ordered<T, E>,
    mut reader: R,
    frames: &Range<usize>,
    block: usize,
    mut decode: F,
) where
    R: TrajectoryReader,
    E: From<TrajectoryError>,
    F: FnMut(&mut R, Range<usize>) -> Result<T, E>,
{
    let _retire = Retire(shared);
    let mut at = frames.start;
    while let Some(id) = shared.claim() {
        let start = frames.start + id * block;
        let end = (start + block).min(frames.end);
        if at != start
            && let Err(error) = reader.seek(start)
        {
            shared.complete(id, Err(error.into()));
            return;
        }
        let result = decode(&mut reader, start..end);
        let failed = result.is_err();
        shared.complete(id, result);
        if failed {
            return;
        }
        at = end;
    }
}

/// Blocks needed to cover `frames` and how far ahead of the consumer they may run.
pub(crate) fn plan(frames: usize, block: usize, workers: usize) -> (usize, usize) {
    (
        frames.div_ceil(block.max(1)),
        workers.saturating_mul(2).max(2),
    )
}

#[cfg(test)]
#[path = "ordered_blocks_tests.rs"]
mod tests;
