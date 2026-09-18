use super::*;
use crate::{MemoryBudget, ScratchPolicy};
use std::sync::atomic::{AtomicUsize, Ordering};

#[test]
fn ordered_consumption_is_bounded_and_identical_across_worker_counts() {
    for workers in [1, 2, 4, 8] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .memory_budget(MemoryBudget::new(2048).expect("budget"))
            .scratch_policy(ScratchPolicy::new(0))
            .build()
            .expect("context");
        let mut output = Vec::new();
        try_for_each_block_in(
            BlockPlan::new(1000, 1),
            &context,
            100,
            |index, _| Ok::<_, std::convert::Infallible>(index),
            |index| {
                output.push(index);
                Ok(())
            },
        )
        .expect("execution");
        assert_eq!(output, (0..1000).collect::<Vec<_>>());
        assert!(context.peak_reserved_bytes() <= 2048);
        assert_eq!(context.reserved_bytes(), 0);
    }
}

#[test]
fn no_block_runs_when_its_working_set_cannot_fit() {
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(64).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let calls = AtomicUsize::new(0);
    let result = try_for_each_block_in(
        BlockPlan::new(1, 1),
        &context,
        128,
        |_, _| {
            calls.fetch_add(1, Ordering::Relaxed);
            Ok::<_, std::convert::Infallible>(())
        },
        |()| Ok(()),
    );
    assert!(matches!(result, Err(BlockExecutionError::Memory(_))));
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn cancellation_in_the_consumer_stops_future_windows_and_releases_memory() {
    let context = ExecutionContext::default();
    let mut count = 0;
    let result = try_for_each_block_in(
        BlockPlan::new(1000, 1),
        &context,
        128,
        |_, _| Ok::<_, std::convert::Infallible>(()),
        |()| {
            count += 1;
            context.cancellation().cancel();
            Ok(())
        },
    );
    assert!(matches!(result, Err(BlockExecutionError::Cancelled)));
    assert_eq!(count, 1);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn nested_blocks_finish_without_a_private_pool_or_deadlock() {
    let context = ExecutionContext::default();
    let mut total = 0;
    try_for_each_block_in(
        BlockPlan::new(32, 1),
        &context,
        32,
        |_, _| {
            let mut inner = 0;
            try_for_each_block_in(
                BlockPlan::new(32, 1),
                &context,
                32,
                |_, _| Ok::<_, std::convert::Infallible>(1),
                |value| {
                    inner += value;
                    Ok(())
                },
            )
            .expect("nested work");
            Ok::<_, std::convert::Infallible>(inner)
        },
        |value| {
            total += value;
            Ok(())
        },
    )
    .expect("outer work");
    assert_eq!(total, 1024);
}

#[test]
fn the_first_result_is_consumed_before_a_later_block_finishes() {
    let context = ExecutionContext::builder()
        .worker_budget(2)
        .build()
        .expect("context");
    if context.worker_budget() < 2 {
        return;
    }
    let (sent, received) = std::sync::mpsc::channel();
    let received = std::sync::Mutex::new(received);
    try_for_each_block_in(
        BlockPlan::new(2, 1),
        &context,
        32,
        |index, _| {
            if index == 1 {
                received
                    .lock()
                    .expect("receiver")
                    .recv_timeout(std::time::Duration::from_secs(5))
                    .expect("first result must be consumed while the second block is running");
            }
            Ok::<_, std::convert::Infallible>(index)
        },
        |index| {
            if index == 0 {
                sent.send(()).expect("consumer notification");
            }
            Ok(())
        },
    )
    .expect("overlapping consumption");
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn worker_panics_release_all_pending_results_and_reservations() {
    let context = ExecutionContext::default();
    let result = try_for_each_block_in(
        BlockPlan::new(32, 1),
        &context,
        128,
        |index, _| {
            assert_ne!(index, 1, "injected worker failure");
            Ok::<_, std::convert::Infallible>(vec![0_u8; 128])
        },
        |_| Ok(()),
    );
    assert!(matches!(result, Err(BlockExecutionError::Worker(_))));
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn a_nested_consumer_executes_inline_when_the_outer_window_owns_all_admission() {
    let context = ExecutionContext::builder()
        .worker_budget(1)
        .build()
        .expect("context");
    let mut total = 0;
    try_for_each_block_in(
        BlockPlan::new(4, 1),
        &context,
        8,
        |_, _| Ok::<_, std::convert::Infallible>(()),
        |()| {
            try_for_each_block_in(
                BlockPlan::new(4, 1),
                &context,
                8,
                |_, _| Ok::<_, std::convert::Infallible>(1),
                |value| {
                    total += value;
                    Ok(())
                },
            )
            .expect("nested consumer");
            Ok(())
        },
    )
    .expect("outer consumer");
    assert_eq!(total, 16);
    assert_eq!(context.reserved_bytes(), 0);
}
