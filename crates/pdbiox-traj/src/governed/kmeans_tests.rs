use super::*;

#[test]
fn borrowed_kmeans_preserves_governed_coverage() {
    let policy = AnalysisPolicy::default();
    let values = [0.0, 0.2, 5.0, 5.2];
    let Ok(result) = analyse_kmeans_view(
        &values,
        4,
        1,
        KMeansOptions {
            initial_centres: &[0, 2],
            maximum_iterations: 20,
            convergence_tolerance_squared: 1e-12,
        },
        &policy,
    ) else {
        panic!("valid governed borrowed k-means");
    };
    assert_eq!(result.coverage.used, 4);
    assert_eq!(result.value.labels, vec![0, 0, 1, 1]);
}
