#![allow(
    clippy::float_cmp,
    reason = "exact identities on exactly-representable inputs"
)]
use super::*;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn a_diagonal_matrix_reports_its_diagonal_largest_first() {
    let decomposition = symmetric([[1.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 3.0]]);
    assert!(close(decomposition.values[0], 5.0));
    assert!(close(decomposition.values[1], 3.0));
    assert!(close(decomposition.values[2], 1.0));
}

#[test]
fn eigenvectors_satisfy_the_relation_that_defines_them() {
    let matrix = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]];
    let decomposition = symmetric(matrix);
    for position in 0..3 {
        let Some(vector) = decomposition.vector(position) else {
            panic!("expected three vectors")
        };
        let value = decomposition.values[position];
        for row in 0..3 {
            let mut product = 0.0;
            for column in 0..3 {
                product += matrix[row][column] * vector[column];
            }
            assert!(
                (product - value * vector[row]).abs() < 1e-8,
                "row {row} of eigenvector {position}"
            );
        }
    }
}

#[test]
fn eigenvectors_are_unit_length() {
    let decomposition = symmetric([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]]);
    for position in 0..3 {
        let Some(vector) = decomposition.vector(position) else {
            panic!("expected a vector")
        };
        let length = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
        assert!(
            (length - 1.0).abs() < 1e-9,
            "vector {position} has length {length}"
        );
    }
}

#[test]
fn the_dominant_vector_is_the_one_belonging_to_the_largest_value() {
    let decomposition = symmetric([[1.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 2.0]]);
    let dominant = decomposition.dominant();
    assert!(
        close(dominant[1].abs(), 1.0),
        "the second axis carries the largest value"
    );
}

#[test]
fn the_same_matrix_always_decomposes_the_same_way() {
    let matrix = [[2.0, -1.0, 0.3], [-1.0, 2.0, 0.1], [0.3, 0.1, 1.0]];
    let first = symmetric(matrix);
    let second = symmetric(matrix);
    assert_eq!(first.values, second.values);
    assert_eq!(first.vectors, second.vectors);
}

#[test]
fn a_four_by_four_decomposes_as_readily_as_a_three_by_three() {
    let decomposition = symmetric([
        [4.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 7.0, 0.0],
        [0.0, 0.0, 0.0, 2.0],
    ]);
    assert!(close(decomposition.values[0], 7.0));
    assert!(close(decomposition.values[3], 1.0));
}

#[test]
fn a_position_past_the_dimension_yields_no_vector() {
    let decomposition = symmetric([[1.0, 0.0], [0.0, 2.0]]);
    assert!(decomposition.vector(2).is_none());
    assert!(decomposition.vector(0).is_some());
}
