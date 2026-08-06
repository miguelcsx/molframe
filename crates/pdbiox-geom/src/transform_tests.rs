#![allow(clippy::float_cmp, reason = "the identity is exact by construction")]

use super::*;

#[test]
fn the_identity_leaves_a_position_where_it_was() {
    let point = [1.5, -2.0, 0.25];
    assert_eq!(Rigid::IDENTITY.apply(point), point);
    assert_eq!(Rigid::default().apply(point), point);
}

#[test]
fn a_transform_undone_returns_the_original_position() {
    let quarter_turn = Rigid::new(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [3.0, -4.0, 5.0],
    );
    let point = [1.0, 2.0, 3.0];
    let there_and_back = quarter_turn.inverse().apply(quarter_turn.apply(point));
    for axis in 0..3 {
        assert!((there_and_back[axis] - point[axis]).abs() < 1e-5);
    }
}

#[test]
fn composing_two_transforms_matches_applying_them_in_turn() {
    let first = Rigid::translation([1.0, 0.0, 0.0]);
    let second = Rigid::new(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [0.0, 0.0, 2.0],
    );
    let point = [1.0, 2.0, 3.0];
    let stepwise = second.apply(first.apply(point));
    let composed = first.then(&second).apply(point);
    for axis in 0..3 {
        assert!(
            (stepwise[axis] - composed[axis]).abs() < 1e-5,
            "axis {axis}"
        );
    }
}

#[test]
fn a_rotation_preserves_distances() {
    let turn = Rigid::new(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [7.0, 8.0, 9.0],
    );
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 0.0, -1.0];
    let before = crate::distance(a, b);
    let after = crate::distance(turn.apply(a), turn.apply(b));
    assert!((before - after).abs() < 1e-5);
}

#[test]
fn a_rotation_has_unit_determinant_and_a_reflection_does_not() {
    assert!((Rigid::IDENTITY.determinant() - 1.0).abs() < 1e-12);
    let mirror = Rigid::new(
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [0.0; 3],
    );
    assert!((mirror.determinant() + 1.0).abs() < 1e-12);
}

#[test]
fn applying_to_a_whole_set_matches_applying_one_at_a_time() {
    let turn = Rigid::new(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [1.0, 1.0, 1.0],
    );
    let original = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut batch = original;
    turn.apply_all(&mut batch);
    for (index, point) in original.iter().enumerate() {
        assert_eq!(batch[index], turn.apply(*point));
    }
}
