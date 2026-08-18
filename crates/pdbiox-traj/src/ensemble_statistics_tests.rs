use super::{
    RemainderPolicy, block_convergence, group_coordinate_variance, group_coordinate_variance_view,
};

#[test]
fn group_variance_uses_explicit_grouping_and_current_frame() {
    let frames = vec![
        vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
        vec![[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
    ];
    let Ok(groups) = group_coordinate_variance(&frames, &[vec![0], vec![1]]) else {
        panic!("valid groups");
    };
    assert_vector_close(groups[0].variance_by_axis, [1.0, 0.0, 0.0]);
    assert_vector_close(groups[1].variance_by_axis, [0.0; 3]);
}

#[test]
fn borrowed_group_variance_matches_owned_frames() {
    let positions = [
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
    ];
    let view = crate::FrameView::new(&positions, 2, 2).expect("valid borrowed frames");
    let owned = vec![
        vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
        vec![[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]],
    ];
    let groups = [vec![0], vec![1]];
    assert_eq!(
        group_coordinate_variance_view(view, &groups),
        group_coordinate_variance(&owned, &groups)
    );
}

#[test]
fn block_remainder_policy_is_not_implicit() {
    let values = [1.0, 3.0, 5.0];
    assert!(block_convergence(&values, 2, RemainderPolicy::Reject).is_err());
    let Ok(blocks) = block_convergence(&values, 2, RemainderPolicy::Include) else {
        panic!("explicit partial block");
    };
    assert!((blocks[0].block_mean - 2.0).abs() < 1.0e-12);
    assert!((blocks[1].block_mean - 5.0).abs() < 1.0e-12);
    assert!((blocks[1].cumulative_mean - 3.0).abs() < 1.0e-12);
}

fn assert_vector_close(actual: [f64; 3], expected: [f64; 3]) {
    actual
        .into_iter()
        .zip(expected)
        .for_each(|(actual, expected)| assert!((actual - expected).abs() < 1.0e-12));
}
