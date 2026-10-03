use super::*;

fn diagonal(values: [f64; 4]) -> [[f64; 4]; 4] {
    let mut matrix = [[0.0; 4]; 4];
    for (index, value) in values.into_iter().enumerate() {
        matrix[index][index] = value;
    }
    matrix
}

#[test]
fn a_diagonal_matrix_yields_its_largest_entry_and_axis() {
    let k = diagonal([-1.0, 3.0, -1.0, -1.0]);
    let Some((lambda, vector)) = dominant_eigenpair(&k, 4.0) else {
        panic!("expected a pair")
    };
    assert!((lambda - 3.0).abs() < 1e-10);
    assert!((vector[1].abs() - 1.0).abs() < 1e-10);
    assert_eq!(
        max_eigenvalue(&k, 4.0).map(|value| (value - 3.0).abs() < 1e-10),
        Some(true)
    );
}

#[test]
fn a_repeated_top_eigenvalue_is_refused() {
    let k = diagonal([1.0, 1.0, 1.0, -3.0]);
    assert!(dominant_eigenpair(&k, 3.0).is_none());
}

#[test]
fn a_non_positive_start_is_refused() {
    assert!(max_eigenvalue(&diagonal([0.0; 4]), 0.0).is_none());
}
