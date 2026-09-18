//! Governed trajectory and ensemble analyses.

mod dimensionality;
mod ensemble;
mod kmeans;

pub use dimensionality::{
    analyse_cartesian_pca, analyse_cartesian_pca_view, analyse_diffusion_map, analyse_dihedral_pca,
};
pub use ensemble::{
    GovernedEnsembleError, analyse_agglomerative_clustering, analyse_block_convergence,
    analyse_cluster_population_similarity, analyse_dbscan_clustering,
    analyse_generalized_procrustes_mean, analyse_generalized_procrustes_mean_view,
    analyse_group_coordinate_variance, analyse_group_coordinate_variance_view,
    analyse_harmonic_ensemble_similarity, analyse_mean_squared_displacement_view, analyse_medoid,
    analyse_pairwise_fitted_rmsd, analyse_pairwise_fitted_rmsd_view,
    analyse_pairwise_torus_distance, analyse_rmsd_to_reference, analyse_rmsd_to_reference_view,
};
pub use kmeans::{analyse_kmeans, analyse_kmeans_view};
