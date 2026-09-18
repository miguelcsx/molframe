use super::*;
use crate::ExecutionContext;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

struct ObserveDrop {
    context: ExecutionContext,
    observed: Arc<AtomicUsize>,
}

impl Drop for ObserveDrop {
    fn drop(&mut self) {
        self.observed
            .store(self.context.reserved_bytes(), Ordering::SeqCst);
    }
}

#[test]
fn rejection_drops_the_value_before_releasing_its_reservation() {
    let context = ExecutionContext::default();
    let observed = Arc::new(AtomicUsize::new(0));
    let value = ObserveDrop {
        context: context.clone(),
        observed: observed.clone(),
    };
    assert!(Retained::new(value, context.try_reserve(8).expect("reservation"), 16).is_err());
    assert_eq!(observed.load(Ordering::SeqCst), 8);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn shared_results_keep_the_charge_until_the_final_owner_is_dropped() {
    let context = ExecutionContext::default();
    let observed = Arc::new(AtomicUsize::new(0));
    let retained = Retained::new(
        ObserveDrop {
            context: context.clone(),
            observed: observed.clone(),
        },
        context.try_reserve(16).expect("reservation"),
        8,
    )
    .expect("retained value");
    let first = Arc::new(retained);
    let second = first.clone();
    drop(first);
    assert_eq!(context.reserved_bytes(), 8);
    drop(second);
    assert_eq!(observed.load(Ordering::SeqCst), 8);
    assert_eq!(context.reserved_bytes(), 0);
}
