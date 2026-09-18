use super::{
    HarmonicSimilarityOptions, cluster_population_similarity, harmonic_ensemble_similarity,
};

#[test]
fn identical_harmonic_ensembles_have_unit_similarity() {
    let observations = vec![vec![-1.0, 0.0], vec![1.0, 0.0]];
    let Ok(similarity) = harmonic_ensemble_similarity(
        &observations,
        &observations,
        HarmonicSimilarityOptions {
            covariance_regularization: 1e-6,
            memory_limit_bytes: 1_024,
        },
    ) else {
        panic!("valid harmonic ensemble");
    };
    assert!((similarity.similarity - 1.0).abs() < 1e-10);
}

#[test]
fn cluster_similarity_spans_identical_and_disjoint_populations() {
    assert_eq!(cluster_population_similarity(&[0, 1], &[0, 1], 2), Ok(1.0));
    let Ok(disjoint) = cluster_population_similarity(&[0, 0], &[1, 1], 2) else {
        panic!("valid populations");
    };
    assert!(disjoint.abs() < 1e-12);
}

#[test]
fn regularization_and_memory_are_never_hidden_defaults() {
    let observations = vec![vec![0.0], vec![1.0]];
    assert!(
        harmonic_ensemble_similarity(
            &observations,
            &observations,
            HarmonicSimilarityOptions {
                covariance_regularization: 0.0,
                memory_limit_bytes: 1_024,
            },
        )
        .is_err()
    );
}
