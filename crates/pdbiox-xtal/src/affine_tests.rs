use super::AffineTransform;

#[test]
fn a_general_matrix_is_preserved_and_applied() {
    let transform = AffineTransform::new(
        [[2.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [1.0, 2.0, 3.0],
    );
    assert_near(transform.apply([1.0, 2.0, 4.0]), [4.0, 4.0, -1.0]);
}

#[test]
fn composition_matches_stepwise_application() {
    let first = AffineTransform::new(
        [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [1.0, 0.0, 0.0],
    );
    let second = AffineTransform::new(AffineTransform::IDENTITY.matrix, [0.0, 3.0, 0.0]);
    let point = [2.0, 1.0, 0.0];
    assert_near(
        first.then(&second).apply(point),
        second.apply(first.apply(point)),
    );
}

fn assert_near(actual: [f32; 3], expected: [f32; 3]) {
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (*actual - expected).abs() < 1e-6)
    );
}
