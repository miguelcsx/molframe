//! Governed trajectory and ensemble analyses.

mod dimensionality;
mod ensemble;

pub use dimensionality::{analyse_cartesian_pca, analyse_diffusion_map, analyse_dihedral_pca};
pub use ensemble::{
    GovernedEnsembleError, analyse_agglomerative_clustering, analyse_block_convergence,
    analyse_cluster_population_similarity, analyse_dbscan_clustering,
    analyse_generalized_procrustes_mean, analyse_group_coordinate_variance,
    analyse_harmonic_ensemble_similarity, analyse_kmeans, analyse_medoid,
    analyse_pairwise_fitted_rmsd, analyse_pairwise_torus_distance, analyse_rmsd_to_reference,
};
