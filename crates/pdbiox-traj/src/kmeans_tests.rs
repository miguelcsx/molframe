use super::{KMeansOptions, kmeans};

#[test]
fn explicit_initial_centres_find_two_feature_groups() {
    let observations = vec![vec![0.0], vec![0.2], vec![5.0], vec![5.2]];
    let Ok(result) = kmeans(
        &observations,
        KMeansOptions {
            initial_centres: &[0, 2],
            maximum_iterations: 20,
            convergence_tolerance_squared: 1e-12,
        },
    ) else {
        panic!("valid k-means");
    };
    assert_eq!(result.labels, vec![0, 0, 1, 1]);
    assert_eq!(result.centres, vec![vec![0.1], vec![5.1]]);
}

#[test]
fn random_or_duplicate_initialization_is_not_inferred() {
    let result = kmeans(
        &[vec![0.0], vec![1.0]],
        KMeansOptions {
            initial_centres: &[0, 0],
            maximum_iterations: 10,
            convergence_tolerance_squared: 0.0,
        },
    );
    assert!(result.is_err());
}
