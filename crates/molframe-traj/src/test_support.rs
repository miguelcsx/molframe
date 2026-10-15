//! Fixtures shared by unit tests across modules.

use molframe_core::{ExecutionContext, MemoryBudget, ScratchPolicy};

/// Execution context with the given worker budget.
pub(crate) fn worker_context(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("context")
}

/// Execution context with a memory budget of `bytes` and no scratch pool.
pub(crate) fn budget_context(bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(bytes).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context")
}

/// Asserts two scalars agree to within single-precision test tolerance.
pub(crate) fn assert_close(actual: f32, expected: f32) {
    assert!((actual - expected).abs() < 1.0e-5);
}

/// Asserts two vectors agree component-wise under [`assert_close`].
pub(crate) fn assert_vector_close(actual: [f32; 3], expected: [f32; 3]) {
    actual
        .into_iter()
        .zip(expected)
        .for_each(|(actual, expected)| assert_close(actual, expected));
}
