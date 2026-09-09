use super::*;
use pdbiox_core::{MemoryBudget, ScratchPolicy, structure::UnitCell};

#[test]
fn construction_scratch_is_released_but_nodes_stay_charged() {
    let context = ExecutionContext::default();
    let positions = [[0.0; 3], [f32::NAN; 3], [1.0; 3]];
    let tree = KdTree::build_in(
        &positions,
        &[0, 1, 2],
        None,
        KdPeriodicOptions::default(),
        &context,
    )
    .expect("tree");
    assert_eq!(context.reserved_bytes(), 2 * size_of::<Node>());
    assert_eq!(
        context.peak_reserved_bytes(),
        3 * (size_of::<Entry>() + size_of::<Node>())
    );
    drop(tree);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn periodic_streaming_emits_a_target_once_even_when_many_images_match() {
    let context = ExecutionContext::default();
    let periodic = PeriodicBox::from_cell(UnitCell {
        lengths: [2.0; 3],
        angles: [90.0; 3],
    })
    .expect("cell");
    let positions = [[0.0; 3], [0.5; 3], [1.5; 3]];
    let tree = KdTree::build_in(
        &positions,
        &[0, 1, 2],
        Some(&periodic),
        KdPeriodicOptions::default(),
        &context,
    )
    .expect("tree");
    let mut count = 0;
    tree.for_each_candidate(&[0], 3.0, &context, |left, right, distance| {
        assert_eq!(left, 0);
        assert!(right == 1 || right == 2);
        assert!((distance - 0.75).abs() < 1e-6);
        count += 1;
    })
    .expect("visit");
    assert_eq!(count, 2);
    assert_eq!(
        context.reserved_bytes(),
        tree.nodes.capacity() * size_of::<Node>()
    );
}

#[test]
fn failed_admission_and_cancelled_visits_emit_nothing() {
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(1).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    assert!(matches!(
        KdTree::build_in(
            &[[0.0; 3]],
            &[0],
            None,
            KdPeriodicOptions::default(),
            &context
        ),
        Err(SpatialError::Memory(_))
    ));
    let context = ExecutionContext::default();
    let positions = [[0.0; 3], [1.0; 3]];
    let tree = KdTree::build_in(
        &positions,
        &[0, 1],
        None,
        KdPeriodicOptions::default(),
        &context,
    )
    .expect("tree");
    context.cancellation().cancel();
    assert!(matches!(
        tree.for_each_candidate(&[0], 2.0, &context, |_, _, _| panic!(
            "cancelled visitor emitted"
        )),
        Err(SpatialError::Cancelled)
    ));
}
