use super::*;

#[test]
fn default_context_uses_bounded_memory_and_no_spill() {
    let context = ExecutionContext::default();
    assert_eq!(context.memory_budget(), MemoryBudget::default());
    assert_eq!(context.temp_storage_policy(), &TempStoragePolicy::Disabled);
    assert!(context.worker_budget() >= 1);
}

#[test]
fn contexts_share_one_native_pool_and_keep_independent_worker_budgets() {
    let first = ExecutionContext::builder()
        .worker_budget(1)
        .build()
        .expect("first context");
    let second = ExecutionContext::builder()
        .worker_budget(2)
        .build()
        .expect("second context");
    assert!(Arc::ptr_eq(&first.executor, &second.executor));
    assert_eq!(first.worker_budget(), 1);
    assert_eq!(second.worker_budget(), 2.min(second.executor.capacity()));
    assert!(matches!(
        ExecutionContext::builder().worker_budget(0).build(),
        Err(ContextError::ZeroWorkers)
    ));
}

#[test]
fn builder_rejects_scratch_larger_than_the_memory_budget() {
    let budget = MemoryBudget::new(10).expect("valid budget");
    let result = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(11))
        .build();
    assert!(matches!(
        result,
        Err(ContextError::ScratchExceedsMemory { .. })
    ));
}

#[test]
fn cloned_contexts_share_memory_accounting_and_cancellation() {
    let budget = MemoryBudget::new(10).expect("valid budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(1))
        .build()
        .expect("valid context");
    let clone = context.clone();
    let reservation = context.try_reserve(8).expect("reservation fits");
    assert!(clone.try_reserve(3).is_err());
    clone.cancellation().cancel();
    assert!(context.cancellation().is_cancelled());
    drop(reservation);
    assert!(clone.try_reserve(10).is_ok());
}
