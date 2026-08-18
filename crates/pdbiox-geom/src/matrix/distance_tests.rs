use super::*;

#[test]
fn a_square_matrix_is_symmetric_with_an_exact_zero_diagonal() {
    let points = [[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 12.0]];
    let Ok(matrix) = distance_matrix(&points) else {
        panic!("small square matrix must fit");
    };
    assert_eq!((matrix.rows(), matrix.columns()), (3, 3));
    assert_eq!(matrix.get(0, 0), Some(0.0));
    assert_eq!(matrix.get(0, 1), Some(5.0));
    assert_eq!(matrix.get(1, 0), matrix.get(0, 1));
    assert_eq!(matrix.get(0, 2), Some(12.0));
    assert_eq!(matrix.get(3, 0), None);
}

#[test]
fn a_rectangular_matrix_keeps_row_major_shape_and_order() {
    let left = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let right = [[0.0, 1.0, 0.0]];
    let Ok(matrix) = distance_matrix_between(&left, &right) else {
        panic!("small rectangular matrix must fit");
    };
    assert_eq!((matrix.rows(), matrix.columns()), (2, 1));
    let expected = [1.0, 2.0f64.sqrt()];
    assert!(
        matrix
            .as_slice()
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (actual - expected).abs() < f64::EPSILON)
    );
}

#[test]
fn values_can_be_moved_without_changing_row_major_order() {
    let left = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let right = [[0.0, 1.0, 0.0]];
    let Ok(matrix) = distance_matrix_between(&left, &right) else {
        panic!("small rectangular matrix must fit");
    };
    assert_eq!(matrix.into_values(), vec![1.0, 2.0f64.sqrt()]);
}
