use super::*;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn decomposition<const N: usize>(matrix: [[f64; N]; N]) -> Decomposition<N> {
    let Ok(decomposition) = symmetric(matrix) else {
        panic!("standard eigensolver profile did not converge")
    };
    decomposition
}

#[test]
fn a_diagonal_matrix_reports_its_diagonal_largest_first() {
    let decomposition = decomposition([[1.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 3.0]]);
    assert!(close(decomposition.values[0], 5.0));
    assert!(close(decomposition.values[1], 3.0));
    assert!(close(decomposition.values[2], 1.0));
}

#[test]
fn eigenvectors_satisfy_the_relation_that_defines_them() {
    let matrix = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]];
    let decomposition = decomposition(matrix);
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
    let decomposition = decomposition([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]]);
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
    let decomposition = decomposition([[1.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 2.0]]);
    let dominant = decomposition.dominant();
    assert!(
        close(dominant[1].abs(), 1.0),
        "the second axis carries the largest value"
    );
}

#[test]
fn the_same_matrix_always_decomposes_the_same_way() {
    let matrix = [[2.0, -1.0, 0.3], [-1.0, 2.0, 0.1], [0.3, 0.1, 1.0]];
    let first = decomposition(matrix);
    let second = decomposition(matrix);
    assert!(
        first
            .values
            .iter()
            .zip(second.values)
            .all(|(left, right)| left.to_bits() == right.to_bits())
    );
    assert!(
        first
            .vectors
            .iter()
            .flatten()
            .zip(second.vectors.iter().flatten())
            .all(|(left, right)| left.to_bits() == right.to_bits())
    );
}

#[test]
fn a_four_by_four_decomposes_as_readily_as_a_three_by_three() {
    let decomposition = decomposition([
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
    let decomposition = decomposition([[1.0, 0.0], [0.0, 2.0]]);
    assert!(decomposition.vector(2).is_none());
    assert!(decomposition.vector(0).is_some());
}

#[test]
fn invalid_and_insufficient_convergence_controls_are_reported() {
    let matrix = [[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]];
    assert_eq!(
        symmetric_with_options(
            matrix,
            EigenOptions {
                relative_tolerance: 0.0,
                maximum_sweeps: 24,
            },
        ),
        Err(EigenError::InvalidOptions)
    );
    assert_eq!(
        symmetric_with_options(
            matrix,
            EigenOptions {
                relative_tolerance: f64::EPSILON,
                maximum_sweeps: 1,
            },
        ),
        Err(EigenError::DidNotConverge)
    );
}

#[test]
fn non_finite_and_overflowing_matrices_are_reported() {
    assert_eq!(
        symmetric([[f64::NAN, 0.0], [0.0, 1.0]]),
        Err(EigenError::NonFiniteMatrix)
    );
    assert_eq!(
        symmetric([[f64::MAX, 0.0], [0.0, f64::MAX]]),
        Err(EigenError::NonFiniteMatrix)
    );
}
