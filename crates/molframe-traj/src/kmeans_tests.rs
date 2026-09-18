use super::{KMeansError, KMeansOptions, kmeans, kmeans_view};

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

#[test]
fn borrowed_row_major_observations_match_owned_kmeans() {
    let owned = vec![vec![0.0], vec![0.2], vec![5.0], vec![5.2]];
    let values = [0.0, 0.2, 5.0, 5.2];
    let options = KMeansOptions {
        initial_centres: &[0, 2],
        maximum_iterations: 20,
        convergence_tolerance_squared: 1e-12,
    };
    let owned_result = kmeans(&owned, options);
    let borrowed_result = kmeans_view(&values, 4, 1, options);
    assert_eq!(borrowed_result, owned_result);
}

#[test]
fn borrowed_observations_require_an_exact_shape() {
    let result = kmeans_view(
        &[0.0, 0.2, 5.0],
        2,
        2,
        KMeansOptions {
            initial_centres: &[0],
            maximum_iterations: 20,
            convergence_tolerance_squared: 1e-12,
        },
    );
    assert_eq!(result, Err(KMeansError::InvalidObservations));
}
