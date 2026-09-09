use super::*;

#[test]
fn default_budget_is_one_hundred_megabytes() {
    assert_eq!(MemoryBudget::default().bytes(), 100_000_000);
}

#[test]
fn budgets_reject_zero_but_accept_any_positive_size() {
    assert_eq!(MemoryBudget::new(0), Err(MemoryBudgetError::Zero));
    assert!(MemoryBudget::new(1).is_ok());
    assert!(
        MemoryBudget::new(64_000_000_000).is_ok(),
        "a caller who has provisioned the machine sets the ceiling, not the library"
    );
    assert!(MemoryBudget::new(usize::MAX).is_ok());
}

#[test]
fn reservations_are_released_when_the_permit_is_dropped() {
    let account = Arc::new(MemoryAccount::new(
        MemoryBudget::new(10).expect("valid budget"),
    ));
    let first = MemoryAccount::reserve(&account, 7).expect("reservation fits");
    assert!(MemoryAccount::reserve(&account, 4).is_err());
    drop(first);
    assert!(MemoryAccount::reserve(&account, 10).is_ok());
}

#[test]
fn concurrent_reservations_never_exceed_the_budget() {
    let account = Arc::new(MemoryAccount::new(
        MemoryBudget::new(1).expect("valid budget"),
    ));
    let first = MemoryAccount::reserve(&account, 1).expect("reservation fits");
    let other = Arc::clone(&account);
    let attempt = std::thread::spawn(move || MemoryAccount::reserve(&other, 1).is_err());
    assert!(attempt.join().expect("thread completes"));
    drop(first);
}
