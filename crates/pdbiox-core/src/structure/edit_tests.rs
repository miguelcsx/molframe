use crate::diagnostic::Code;
use crate::index::ModelIndex;
use crate::{ExecutionContext, MemoryBudget, ScratchPolicy};

fn context() -> ExecutionContext {
    ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context")
}

#[test]
fn committing_coordinates_creates_a_new_generation_and_keeps_the_original() {
    let original = crate::structure::fixture::sample();
    let before = original.positions()[0];
    let mut editor = original.edit_coordinates(&context()).expect("edit fits");
    let Some(positions) = editor.positions_mut(ModelIndex::new(0)) else {
        panic!("first model missing")
    };
    positions[0] = [9.0, 8.0, 7.0];
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };

    assert!(
        original.positions()[0]
            .iter()
            .zip(before)
            .all(|(actual, expected)| (actual - expected).abs() < f32::EPSILON)
    );
    assert!(
        edited.positions()[0]
            .iter()
            .zip([9.0, 8.0, 7.0])
            .all(|(actual, expected)| (actual - expected).abs() < f32::EPSILON)
    );
    assert_eq!(original.generation().get(), 0);
    assert_eq!(edited.generation().get(), 1);
}

#[test]
fn a_non_finite_coordinate_aborts_the_transaction() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit_coordinates(&context()).expect("edit fits");
    let Some(positions) = editor.positions_mut(ModelIndex::new(0)) else {
        panic!("first model missing")
    };
    positions[0][0] = f32::NAN;
    let Err(findings) = editor.commit() else {
        panic!("non-finite coordinate committed")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E3008));
    assert!(original.positions()[0][0].is_finite());
}

#[test]
fn checked_coordinate_access_uses_the_diagnostic_registry() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit_coordinates(&context()).expect("edit fits");
    let result = editor.try_positions_mut(ModelIndex::new(1));
    assert_eq!(
        result.err().map(|finding| finding.code()),
        Some(Code::E6003)
    );
}

#[test]
fn budgeted_edit_charges_the_detached_backing_until_every_snapshot_owner_drops() {
    let original = crate::structure::fixture::sample();
    let bytes = original
        .data()
        .coords
        .block(ModelIndex::new(0))
        .expect("model")
        .allocated_bytes();
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(bytes).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let editor = original.edit_coordinates(&context).expect("copy fits");
    assert_eq!(context.reserved_bytes(), bytes);
    let snapshot = editor.snapshot().expect("valid snapshot");
    drop(editor);
    assert_eq!(context.reserved_bytes(), bytes);
    let clone = snapshot.clone();
    drop(snapshot);
    assert_eq!(context.reserved_bytes(), bytes);
    drop(clone);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn budgeted_edit_refuses_before_detaching_coordinates() {
    let original = crate::structure::fixture::sample();
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(1).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    assert!(original.edit_coordinates(&context).is_err());
    assert_eq!(context.reserved_bytes(), 0);
}
