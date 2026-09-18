use super::*;
use crate::{MemoryBudget, ScratchPolicy};
use std::sync::Arc;

#[derive(Debug)]
struct TestBatch {
    rows: usize,
    bytes: usize,
}

impl Batch for TestBatch {
    fn rows(&self) -> usize {
        self.rows
    }

    fn retained_bytes(&self) -> usize {
        self.bytes
    }
}

#[test]
fn zero_demand_applies_backpressure() {
    assert!(!BatchDemand::new(0, 1).can_accept_work());
    assert!(!BatchDemand::new(1, 0).can_accept_work());
}

#[test]
fn demand_checks_both_rows_and_retained_bytes() {
    let demand = BatchDemand::new(4, 16);
    assert!(demand.accepts(&TestBatch { rows: 4, bytes: 16 }));
    assert!(!demand.accepts(&TestBatch { rows: 5, bytes: 16 }));
    assert!(!demand.accepts(&TestBatch { rows: 4, bytes: 17 }));
}

#[test]
fn a_batch_lease_holds_capacity_until_it_is_dropped() {
    let budget = MemoryBudget::new(16).expect("valid budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(1))
        .build()
        .expect("valid context");
    let lease =
        BatchLease::try_new(TestBatch { rows: 4, bytes: 12 }, &context).expect("batch fits");
    assert_eq!(lease.batch().rows(), 4);
    assert_eq!(context.reserved_bytes(), 12);
    assert_eq!(context.peak_reserved_bytes(), 12);
    assert_eq!(context.live_batches(), 1);
    assert_eq!(context.peak_live_batches(), 1);
    assert!(context.try_reserve(5).is_err());
    drop(lease);
    assert_eq!(context.reserved_bytes(), 0);
    assert_eq!(context.live_batches(), 0);
    assert!(context.try_reserve(16).is_ok());
}

#[test]
fn precharged_lease_releases_unused_slot_capacity() {
    let budget = MemoryBudget::new(32).expect("valid budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context");
    let reservation = context.try_reserve(24).expect("slot fits");
    let lease = BatchLease::try_from_reservation(TestBatch { rows: 1, bytes: 8 }, reservation)
        .expect("batch fits slot");
    assert!(context.try_reserve(24).is_ok());
    drop(lease);
    assert!(context.try_reserve(32).is_ok());
}

#[test]
fn recycling_keeps_capacity_charged_and_ends_the_live_batch() {
    let budget = MemoryBudget::new(32).expect("valid budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context");
    let pool = Arc::new(BatchBufferPool::new());
    let reservation = context.try_reserve(24).expect("slot fits");
    let lease = BatchLease::try_from_reservation_recycling(
        TestBatch { rows: 1, bytes: 8 },
        reservation,
        Arc::clone(&pool),
    )
    .expect("batch fits slot");
    assert_eq!(context.reserved_bytes(), 8);
    assert_eq!(context.live_batches(), 1);
    drop(lease);
    assert_eq!(context.reserved_bytes(), 8);
    assert_eq!(context.live_batches(), 0);
    let (batch, reservation) = pool.take().expect("recycled batch");
    assert_eq!(batch.rows(), 1);
    assert_eq!(reservation.bytes(), 8);
    drop(reservation);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn recycling_retains_only_one_of_concurrent_batches() {
    let budget = MemoryBudget::new(32).expect("valid budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("valid context");
    let pool = Arc::new(BatchBufferPool::new());
    let first = BatchLease::try_from_reservation_recycling(
        TestBatch { rows: 1, bytes: 8 },
        context.try_reserve(8).expect("first slot"),
        Arc::clone(&pool),
    )
    .expect("first lease");
    let second = BatchLease::try_from_reservation_recycling(
        TestBatch { rows: 1, bytes: 8 },
        context.try_reserve(8).expect("second slot"),
        Arc::clone(&pool),
    )
    .expect("second lease");
    drop(first);
    drop(second);
    assert_eq!(context.reserved_bytes(), 8);
    drop(pool);
    assert_eq!(context.reserved_bytes(), 0);
}

struct ObservedBatch {
    context: ExecutionContext,
    observed: Arc<std::sync::atomic::AtomicUsize>,
    bytes: usize,
}

impl Batch for ObservedBatch {
    fn rows(&self) -> usize {
        1
    }
    fn retained_bytes(&self) -> usize {
        self.bytes
    }
}

impl Drop for ObservedBatch {
    fn drop(&mut self) {
        self.observed.store(
            self.context.reserved_bytes(),
            std::sync::atomic::Ordering::SeqCst,
        );
    }
}

#[test]
fn rejected_recycling_keeps_the_charge_until_the_buffer_is_destroyed() {
    let context = ExecutionContext::default();
    let observed = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let pool = Arc::new(BatchBufferPool::new());
    for _ in 0..2 {
        let lease = BatchLease::try_from_reservation_recycling(
            ObservedBatch {
                context: context.clone(),
                observed: observed.clone(),
                bytes: 8,
            },
            context.try_reserve(8).expect("reservation"),
            pool.clone(),
        )
        .expect("lease");
        drop(lease);
    }
    assert_eq!(observed.load(std::sync::atomic::Ordering::SeqCst), 16);
    assert_eq!(context.reserved_bytes(), 8);
    drop(pool);
    assert_eq!(observed.load(std::sync::atomic::Ordering::SeqCst), 8);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn an_undersized_reservation_remains_live_during_rejected_batch_destruction() {
    let context = ExecutionContext::default();
    let observed = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let result = BatchLease::try_from_reservation(
        ObservedBatch {
            context: context.clone(),
            observed: observed.clone(),
            bytes: 16,
        },
        context.try_reserve(8).expect("reservation"),
    );
    assert!(result.is_err());
    assert_eq!(observed.load(std::sync::atomic::Ordering::SeqCst), 8);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn sharing_the_lease_keeps_both_capacity_and_accounting_alive() {
    let context = ExecutionContext::default();
    let pool = Arc::new(BatchBufferPool::new());
    let lease = Arc::new(
        BatchLease::try_from_reservation_recycling(
            TestBatch { rows: 1, bytes: 8 },
            context.try_reserve(8).expect("reservation"),
            pool.clone(),
        )
        .expect("lease"),
    );
    let retained = lease.clone();
    drop(lease);
    assert!(pool.take().is_none());
    drop(pool);
    assert_eq!(context.reserved_bytes(), 8);
    assert_eq!(retained.batch().rows(), 1);
    drop(retained);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn the_pool_never_returns_a_batch_while_a_lease_is_still_releasing_it() {
    let context = ExecutionContext::default();
    let pool = BatchBufferPool::new();
    let batch = Arc::new(TestBatch { rows: 1, bytes: 8 });
    assert!(
        pool.put_or_return(batch.clone(), context.try_reserve(8).expect("reservation"))
            .is_none()
    );
    assert!(pool.take().is_none());
    assert_eq!(context.reserved_bytes(), 8);
    drop(batch);
    let (batch, reservation) = pool.take().expect("exclusive batch");
    assert_eq!(batch.rows(), 1);
    drop(reservation);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn the_pool_rejects_an_undercharged_buffer_before_retaining_it() {
    let context = ExecutionContext::default();
    let observed = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let pool = BatchBufferPool::new();
    let result = pool.put(
        ObservedBatch {
            context: context.clone(),
            observed: observed.clone(),
            bytes: 16,
        },
        context.try_reserve(8).expect("reservation"),
    );
    assert!(result.is_err());
    assert!(pool.take().is_none());
    assert_eq!(observed.load(std::sync::atomic::Ordering::SeqCst), 8);
    assert_eq!(context.reserved_bytes(), 0);
}
