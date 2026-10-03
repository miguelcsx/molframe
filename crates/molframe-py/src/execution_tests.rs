use super::*;
use std::time::Instant;

#[test]
fn a_finished_operation_returns_without_waiting_for_a_poll() {
    Python::initialize();
    Python::attach(|py| {
        let started = Instant::now();
        let value = run(py, None, |_| 7).expect("the operation runs");
        assert_eq!(value, 7);
        assert!(
            started.elapsed() < POLL,
            "a quick call must not sleep a whole poll"
        );
    });
}

#[test]
fn cancelling_the_context_stops_an_operation_that_checks_its_token() {
    Python::initialize();
    Python::attach(|py| {
        let context = Bound::new(
            py,
            PyExecutionContext::new(Some(2), None, 0, None, 0).expect("valid limits"),
        )
        .expect("binds");
        let handle = context.clone().unbind();
        std::thread::scope(|scope| {
            scope.spawn(|| {
                std::thread::sleep(Duration::from_millis(120));
                Python::attach(|py| handle.bind(py).get().cancel());
            });
            let stopped = run(py, Some(&*context.borrow()), |context| {
                while !context.cancellation().is_cancelled() {
                    std::thread::sleep(Duration::from_millis(5));
                }
                true
            })
            .expect("a cancelled operation still returns its own result");
            assert!(stopped);
        });
        assert!(context.get().is_cancelled());
    });
}

#[test]
fn limits_the_kernel_cannot_honour_are_refused_when_the_context_is_made() {
    Python::initialize();
    assert!(PyExecutionContext::new(Some(0), None, 0, None, 0).is_err());
    assert!(PyExecutionContext::new(None, Some(0), 0, None, 0).is_err());
    assert!(PyExecutionContext::new(None, Some(10), 1_000, None, 0).is_err());
}

#[test]
fn the_built_context_carries_the_limits_it_was_given() {
    Python::initialize();
    let context = PyExecutionContext::new(Some(3), Some(1 << 20), 0, None, 0).expect("valid");
    let built = context.build(&CancellationToken::new()).expect("builds");
    assert_eq!(built.worker_budget(), 3);
    assert_eq!(built.memory_budget().bytes(), 1 << 20);
}
