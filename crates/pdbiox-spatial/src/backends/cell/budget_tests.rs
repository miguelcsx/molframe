use super::*;
use pdbiox_core::{MemoryBudget, ScratchPolicy};

#[test]
fn an_index_keeps_its_charge_until_its_storage_is_dropped() {
    let context = ExecutionContext::default();
    let positions = [[0.0; 3], [1.0; 3], [2.0; 3]];
    let options = CellGridOptions::default();
    let required = workspace_bytes(&positions, &[0, 1, 2], 1.0, None, options).expect("estimate");
    let index =
        CellList::build_in(&positions, &[0, 1, 2], 1.0, None, options, &context).expect("index");
    assert_eq!(context.peak_reserved_bytes(), required);
    assert!(context.reserved_bytes() < required);
    drop(index);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn a_tiny_budget_refuses_the_grid_before_it_allocates_storage() {
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(1).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let result = CellList::build_in(
        &[[0.0; 3]],
        &[0],
        1.0,
        None,
        CellGridOptions::default(),
        &context,
    );
    assert!(matches!(result, Err(SpatialError::Memory(_))));
    assert_eq!(context.reserved_bytes(), 0);
}
