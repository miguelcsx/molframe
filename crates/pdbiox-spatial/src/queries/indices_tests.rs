use super::*;
use pdbiox_core::{MemoryBudget, ScratchPolicy};

fn context(bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(bytes).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context")
}

#[test]
fn identical_ranges_expand_once_and_keep_the_charge_until_drop() {
    let context = context(16);
    let left = AtomSelection::range(0..4);
    let right = left.clone();
    let indices = IndexWorkspace::new(&left, &right, 4, &context).expect("workspace");
    assert_eq!(indices.left().as_ptr(), indices.right().as_ptr());
    assert_eq!(context.reserved_bytes(), 16);
    drop(indices);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn different_expansions_share_one_budget_and_reject_before_allocation() {
    let context = context(16);
    let left = AtomSelection::range(0..4);
    let right = AtomSelection::range(1..4);
    let result = IndexWorkspace::new(&left, &right, 4, &context);
    assert!(matches!(result, Err(SpatialError::Memory(_))));
    assert_eq!(context.peak_reserved_bytes(), 0);
}

#[test]
fn sparse_selections_borrow_without_charging_another_copy() {
    let context = context(1);
    let left = AtomSelection::Sparse(vec![0, 2, 4]);
    let right = AtomSelection::Sparse(vec![1, 3]);
    let workspace = IndexWorkspace::new(&left, &right, 5, &context).expect("workspace");
    assert_eq!(workspace.left(), &[0, 2, 4]);
    assert_eq!(workspace.right(), &[1, 3]);
    assert_eq!(context.peak_reserved_bytes(), 0);
}

#[test]
fn cancellation_and_invalid_indices_release_every_reservation() {
    let context = context(64);
    assert!(matches!(
        IndexWorkspace::new(
            &AtomSelection::range(0..4),
            &AtomSelection::range(0..8),
            4,
            &context
        ),
        Err(SpatialError::AtomOutOfBounds(_))
    ));
    assert_eq!(context.reserved_bytes(), 0);
    context.cancellation().cancel();
    assert!(matches!(
        IndexWorkspace::new(&AtomSelection::All(4), &AtomSelection::All(4), 4, &context),
        Err(SpatialError::Cancelled)
    ));
    assert_eq!(context.reserved_bytes(), 0);
}
