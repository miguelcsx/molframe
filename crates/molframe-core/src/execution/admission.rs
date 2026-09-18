//! Shared admission for active and pending blocks in one execution graph.

use super::CancellationToken;
use std::cell::Cell;
use std::sync::{Condvar, Mutex};
use std::time::Duration;

thread_local! {
    static DEPTH: Cell<usize> = const { Cell::new(0) };
}

pub(crate) struct ExecutionDepth;

impl ExecutionDepth {
    pub(crate) fn enter() -> (Self, bool) {
        let nested = DEPTH.with(|depth| {
            let previous = depth.get();
            depth.set(previous + 1);
            previous != 0
        });
        (Self, nested || rayon::current_thread_index().is_some())
    }
}

impl Drop for ExecutionDepth {
    fn drop(&mut self) {
        DEPTH.with(|depth| depth.set(depth.get() - 1));
    }
}

#[derive(Debug)]
pub(crate) struct Admission {
    used: Mutex<usize>,
    changed: Condvar,
    limit: usize,
}

impl Admission {
    pub(crate) const fn new(limit: usize) -> Self {
        Self {
            used: Mutex::new(0),
            changed: Condvar::new(),
            limit,
        }
    }

    // Only outer callers wait here, before reserving window storage. Nested
    // work executes inline under its parent's admission instead.
    pub(crate) fn acquire(
        &self,
        requested: usize,
        cancellation: &CancellationToken,
    ) -> Option<AdmissionGuard<'_>> {
        let mut used = match self.used.lock() {
            Ok(used) => used,
            Err(poisoned) => poisoned.into_inner(),
        };
        loop {
            if cancellation.is_cancelled() {
                return None;
            }
            let count = requested.min(self.limit - *used);
            if count > 0 {
                *used += count;
                return Some(AdmissionGuard {
                    admission: self,
                    count,
                });
            }
            let waited = self.changed.wait_timeout(used, Duration::from_millis(50));
            used = match waited {
                Ok((used, _)) => used,
                Err(poisoned) => poisoned.into_inner().0,
            };
        }
    }
}

pub(crate) struct AdmissionGuard<'a> {
    admission: &'a Admission,
    pub(crate) count: usize,
}

impl Drop for AdmissionGuard<'_> {
    fn drop(&mut self) {
        let mut used = match self.admission.used.lock() {
            Ok(used) => used,
            Err(poisoned) => poisoned.into_inner(),
        };
        *used -= self.count;
        self.admission.changed.notify_all();
    }
}

#[cfg(test)]
#[path = "admission_tests.rs"]
mod tests;
