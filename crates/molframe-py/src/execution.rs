//! The execution context Python hands to heavy operations.
//!
//! There is one concept of resource governance in `MolFrame` — the Rust
//! `ExecutionContext` — and this is its Python face: a worker count, a memory
//! budget, scratch and temporary-storage limits, and a cancellation switch.
//! Every heavy binding takes `context=` and runs through [`run`], which also
//! makes it interruptible: Ctrl-C reaches the kernel as cancellation and comes
//! back as `KeyboardInterrupt`.

use molframe::{
    CancellationToken, ExecutionContext, MemoryBudget, ScratchPolicy, TempStoragePolicy,
};
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

/// How often the waiting thread looks for a signal or a cancelled context.
const POLL: Duration = Duration::from_millis(50);

/// Limits and a cancellation switch shared by the operations run under it.
///
/// A context describes limits; each operation builds its own accounting from
/// them, so one context can govern many calls, and `cancel()` stops whichever
/// are running.
#[derive(Debug)]
#[pyclass(
    name = "ExecutionContext",
    frozen,
    skip_from_py_object,
    module = "molframe"
)]
pub struct PyExecutionContext {
    workers: Option<usize>,
    memory_budget: Option<usize>,
    scratch_bytes: usize,
    temp_directory: Option<PathBuf>,
    temp_bytes: u64,
    image_limit: Option<usize>,
    token: CancellationToken,
}

#[pymethods]
impl PyExecutionContext {
    /// Creates a context; omitted limits keep the library's bounded defaults.
    #[new]
    #[pyo3(signature = (*, workers=None, memory_budget=None, scratch_bytes=0, temp_directory=None, temp_bytes=0, image_limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        workers: Option<usize>,
        memory_budget: Option<usize>,
        scratch_bytes: usize,
        temp_directory: Option<PathBuf>,
        temp_bytes: u64,
        image_limit: Option<usize>,
    ) -> PyResult<Self> {
        let context = Self {
            workers,
            memory_budget,
            scratch_bytes,
            temp_directory,
            temp_bytes,
            image_limit,
            token: CancellationToken::new(),
        };
        // Refused here, at the construction the caller wrote, not at the first
        // operation that happens to use it.
        context.build(&CancellationToken::new())?;
        Ok(context)
    }

    /// Asks every operation running under this context to stop.
    fn cancel(&self) {
        self.token.cancel();
    }

    /// Whether `cancel()` has been called.
    #[getter]
    fn is_cancelled(&self) -> bool {
        self.token.is_cancelled()
    }

    /// The worker limit, or `None` for the library default.
    #[getter]
    fn workers(&self) -> Option<usize> {
        self.workers
    }

    /// The memory budget in bytes, or `None` for the library default.
    #[getter]
    fn memory_budget(&self) -> Option<usize> {
        self.memory_budget
    }

    /// The most candidate symmetry images a crystal search may examine, or `None` for the default.
    #[getter]
    const fn image_limit(&self) -> Option<usize> {
        self.image_limit
    }

    /// Bytes of reusable scratch the operations may retain.
    #[getter]
    fn scratch_bytes(&self) -> usize {
        self.scratch_bytes
    }

    /// The directory temporary files may use, or `None` for none.
    #[getter]
    fn temp_directory(&self) -> Option<PathBuf> {
        self.temp_directory.clone()
    }

    /// The most temporary-file bytes the operations may write.
    #[getter]
    fn temp_bytes(&self) -> u64 {
        self.temp_bytes
    }

    fn __repr__(&self) -> String {
        format!(
            "ExecutionContext(workers={:?}, memory_budget={:?}, scratch_bytes={}, temp_bytes={}, cancelled={})",
            self.workers,
            self.memory_budget,
            self.scratch_bytes,
            self.temp_bytes,
            self.token.is_cancelled()
        )
    }
}

impl PyExecutionContext {
    /// The switch `cancel()` flips.
    pub(crate) const fn token(&self) -> &CancellationToken {
        &self.token
    }

    /// The Rust context for one operation, cancelled by `token`.
    pub(crate) fn build(&self, token: &CancellationToken) -> PyResult<ExecutionContext> {
        let mut builder = ExecutionContext::builder()
            .scratch_policy(ScratchPolicy::new(self.scratch_bytes))
            .cancellation(token.clone());
        if let Some(workers) = self.workers {
            builder = builder.worker_budget(workers);
        }
        if let Some(images) = self.image_limit {
            builder = builder.image_search_limit(images);
        }
        if let Some(bytes) = self.memory_budget {
            builder =
                builder.memory_budget(MemoryBudget::new(bytes).map_err(crate::error::kernel)?);
        }
        if let Some(directory) = &self.temp_directory {
            builder = builder.temp_storage_policy(TempStoragePolicy::directory(
                directory.clone(),
                self.temp_bytes,
            ));
        }
        builder.build().map_err(crate::error::kernel)
    }
}

/// Runs `work` under `context`, interruptibly, with the GIL released.
///
/// The calling thread — the one Python delivers signals to — only waits and
/// polls; the kernel runs on a scoped thread. A Ctrl-C, or `cancel()` on the
/// context, cancels the operation's own token, the kernel stops at its next
/// check, and a signal comes back as the `KeyboardInterrupt` it was.
pub(crate) fn run<T, F>(
    py: Python<'_>,
    context: Option<&PyExecutionContext>,
    work: F,
) -> PyResult<T>
where
    T: Send,
    F: FnOnce(&ExecutionContext) -> T + Send,
{
    let token = CancellationToken::new();
    let built = match context {
        Some(context) => context.build(&token)?,
        None => ExecutionContext::builder()
            .cancellation(token.clone())
            .build()
            .map_err(crate::error::kernel)?,
    };
    let user = context.map(|context| context.token.clone());
    // A context cancelled before the call refuses it, without waiting for the
    // first poll to notice.
    if user.as_ref().is_some_and(CancellationToken::is_cancelled) {
        token.cancel();
    }
    let finished = AtomicBool::new(false);
    let waiter = std::thread::current();
    let mut interrupted = None;
    let output = std::thread::scope(|scope| {
        let handle = scope.spawn(|| {
            let output = work(&built);
            finished.store(true, Ordering::Release);
            waiter.unpark();
            output
        });
        while !finished.load(Ordering::Acquire) {
            py.detach(|| std::thread::park_timeout(POLL));
            if finished.load(Ordering::Acquire) {
                break;
            }
            if user.as_ref().is_some_and(CancellationToken::is_cancelled) {
                token.cancel();
            }
            if let Err(signal) = py.check_signals() {
                token.cancel();
                interrupted = Some(signal);
                break;
            }
        }
        py.detach(|| handle.join())
    });
    let output = output.map_err(|_| crate::error::internal("a worker thread panicked"))?;
    match interrupted {
        Some(signal) => Err(signal),
        None => Ok(output),
    }
}

#[cfg(test)]
#[path = "execution_tests.rs"]
mod tests;
