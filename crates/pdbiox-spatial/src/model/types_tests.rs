use super::*;

#[test]
fn resource_and_cancellation_failures_do_not_masquerade_as_bad_selection_values() {
    let memory = SpatialError::Memory(pdbiox_core::MemoryBudgetError::Exhausted {
        requested: 100,
        available: 10,
    });
    assert_eq!(memory.into_diagnostic().code(), Code::E1901);
    assert_eq!(
        SpatialError::Cancelled.into_diagnostic().code(),
        Code::E1904
    );
}
