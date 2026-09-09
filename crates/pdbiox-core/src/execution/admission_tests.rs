use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};

#[test]
fn simultaneous_callers_share_one_admission_ceiling() {
    let admission = Admission::new(4);
    let cancellation = CancellationToken::default();
    let active = AtomicUsize::new(0);
    let peak = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..8 {
            scope.spawn(|| {
                for _ in 0..10 {
                    let guard = admission.acquire(3, &cancellation).expect("admission");
                    let current = active.fetch_add(guard.count, Ordering::SeqCst) + guard.count;
                    peak.fetch_max(current, Ordering::SeqCst);
                    std::thread::yield_now();
                    active.fetch_sub(guard.count, Ordering::SeqCst);
                    drop(guard);
                }
            });
        }
    });
    assert!(peak.load(Ordering::SeqCst) <= 4);
    assert_eq!(*admission.used.lock().expect("counter"), 0);
}

#[test]
fn cancellation_interrupts_a_caller_waiting_for_admission() {
    let admission = Admission::new(1);
    let cancellation = CancellationToken::default();
    let _guard = admission
        .acquire(1, &cancellation)
        .expect("first admission");
    std::thread::scope(|scope| {
        let waiting = scope.spawn(|| admission.acquire(1, &cancellation).is_none());
        cancellation.cancel();
        assert!(waiting.join().expect("waiting caller"));
    });
}

#[test]
fn unwinding_releases_admission_and_nesting_state() {
    let admission = Admission::new(2);
    let cancellation = CancellationToken::default();
    let failed = std::panic::catch_unwind(|| {
        let (_depth, nested) = ExecutionDepth::enter();
        assert!(!nested);
        let _guard = admission.acquire(2, &cancellation).expect("admission");
        panic!("injected failure");
    });
    assert!(failed.is_err());
    let (_depth, nested) = ExecutionDepth::enter();
    assert!(!nested);
    assert_eq!(
        admission.acquire(2, &cancellation).expect("released").count,
        2
    );
}
