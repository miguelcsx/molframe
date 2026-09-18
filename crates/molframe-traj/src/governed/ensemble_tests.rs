use super::*;

fn matrix() -> EnsembleDistanceMatrix {
    EnsembleDistanceMatrix {
        size: 3,
        values: vec![0.0, 0.1, 2.0, 0.1, 0.0, 1.9, 2.0, 1.9, 0.0].into_boxed_slice(),
    }
}

#[test]
fn clustering_medoid_and_convergence_are_governed() {
    let policy = AnalysisPolicy::default();
    let Ok(clusters) = analyse_agglomerative_clustering(&matrix(), 2, Linkage::Average, &policy)
    else {
        panic!("valid governed clustering");
    };
    let Ok(representative) = analyse_medoid(&matrix(), &[0, 1], &policy) else {
        panic!("valid governed medoid");
    };
    let Ok(blocks) = analyse_block_convergence(&[1.0, 2.0], 1, RemainderPolicy::Reject, &policy)
    else {
        panic!("valid governed blocks");
    };
    assert_eq!(clusters.coverage.used, 3);
    assert_eq!(representative.value, 0);
    assert_eq!(blocks.value.len(), 2);
}

#[test]
fn both_similarity_families_are_governed() {
    let policy = AnalysisPolicy::default();
    let observations = vec![vec![0.0], vec![1.0]];
    let Ok(harmonic) = analyse_harmonic_ensemble_similarity(
        &observations,
        &observations,
        HarmonicSimilarityOptions {
            covariance_regularization: 1e-6,
            memory_limit_bytes: 1024,
        },
        &policy,
    ) else {
        panic!("valid governed harmonic similarity");
    };
    let Ok(population) = analyse_cluster_population_similarity(&[0, 1], &[0, 1], 2, &policy) else {
        panic!("valid governed population similarity");
    };
    assert!((harmonic.value.similarity - 1.0).abs() < 1e-12);
    assert!((population.value - 1.0).abs() < 1e-12);
}

#[test]
fn borrowed_frame_views_preserve_governed_results_and_provenance() {
    let positions = [
        [0.0_f32, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.1, 0.0, 0.0],
        [1.1, 0.0, 0.0],
        [0.1, 1.0, 0.0],
        [0.2, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.2, 1.0, 0.0],
    ];
    let view = crate::FrameView::new(&positions, 3, 3).expect("contiguous frames");
    let policy = AnalysisPolicy::default();
    let rmsd = analyse_rmsd_to_reference_view(view, 0, FrameAlignment::Rigid, &policy)
        .expect("governed borrowed RMSD");
    let pairwise = analyse_pairwise_fitted_rmsd_view(view, 1024, &policy)
        .expect("governed borrowed pairwise RMSD");
    let mean = analyse_generalized_procrustes_mean_view(view, 1.0e-6, 20, &policy)
        .expect("governed borrowed mean");
    let msd = analyse_mean_squared_displacement_view(view, &[], 2, &policy)
        .expect("governed borrowed MSD");
    let groups = analyse_group_coordinate_variance_view(view, &[vec![0, 1, 2]], &policy)
        .expect("governed borrowed variance");
    assert_eq!(rmsd.coverage.used, 3);
    assert_eq!(pairwise.value.size, 3);
    assert_eq!(mean.value.len(), 3);
    assert_eq!(msd.value.len(), 3);
    assert_eq!(groups.value.len(), 1);
    let Some(algorithm) = msd.provenance.algorithm.as_ref() else {
        panic!("MSD provenance has an algorithm");
    };
    assert_eq!(algorithm.name(), "mean-squared-displacement");
}
