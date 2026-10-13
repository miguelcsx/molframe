use super::*;

#[test]
fn default_budget_scales_with_the_machine_within_its_floor_and_ceiling() {
    let bytes = MemoryBudget::default().bytes();
    assert!(bytes >= DEFAULT_MEMORY_BUDGET_BYTES, "{bytes}");
    assert!(bytes <= MAX_AUTOMATIC_BUDGET_BYTES, "{bytes}");
}

#[test]
fn physical_memory_maps_to_forty_percent_clamped() {
    let gib = 1_usize << 30;
    assert_eq!(
        MemoryBudget::for_physical_memory(0).bytes(),
        DEFAULT_MEMORY_BUDGET_BYTES
    );
    assert_eq!(
        MemoryBudget::for_physical_memory(10 * gib).bytes(),
        10 * gib / 5 * 2
    );
    assert_eq!(
        MemoryBudget::for_physical_memory(usize::MAX).bytes(),
        MAX_AUTOMATIC_BUDGET_BYTES
    );
}

#[test]
fn budget_text_parses_units_and_rejects_nonsense() {
    assert_eq!(MemoryBudget::parse("1000").unwrap().bytes(), 1000);
    assert_eq!(MemoryBudget::parse("2KB").unwrap().bytes(), 2000);
    assert_eq!(MemoryBudget::parse("2MiB").unwrap().bytes(), 2 << 20);
    assert_eq!(MemoryBudget::parse("1.5G").unwrap().bytes(), 3 << 29);
    assert_eq!(MemoryBudget::parse("0"), Err(MemoryBudgetError::Zero));
    assert!(matches!(
        MemoryBudget::parse("lots"),
        Err(MemoryBudgetError::Unparsable)
    ));
}

#[test]
fn exhaustion_names_figures_and_every_knob() {
    let message = MemoryBudgetError::Exhausted {
        requested: 7,
        available: 3,
    }
    .to_string();
    for needle in [
        "7 bytes",
        "3 bytes",
        "memory_budget",
        "--memory-budget",
        "MOLFRAME_MEMORY_BUDGET",
    ] {
        assert!(message.contains(needle), "{needle}: {message}");
    }
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
