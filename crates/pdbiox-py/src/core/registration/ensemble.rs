//! Registration for governed ensemble projections.

use crate::intrinsic::{
    PyClustering, PyConvergenceBlock, PyEnsembleDistanceMatrix, PyFrameAlignment, PyGroupVariance,
    PyHarmonicSimilarity, PyHarmonicSimilarityOptions, PyKMeans, PyKMeansOptions, PyLinkage,
    PyRemainderPolicy, analyse_agglomerative_clustering, analyse_block_convergence,
    analyse_cluster_population_similarity, analyse_dbscan_clustering,
    analyse_generalized_procrustes_mean, analyse_group_coordinate_variance,
    analyse_harmonic_ensemble_similarity, analyse_kmeans, analyse_medoid,
    analyse_pairwise_fitted_rmsd, analyse_pairwise_torus_distance, analyse_rmsd_to_reference,
};
use pyo3::prelude::*;

pub(super) fn register_ensemble(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFrameAlignment>()?;
    module.add_class::<PyLinkage>()?;
    module.add_class::<PyRemainderPolicy>()?;
    module.add_class::<PyEnsembleDistanceMatrix>()?;
    module.add_class::<PyClustering>()?;
    module.add_class::<PyKMeansOptions>()?;
    module.add_class::<PyKMeans>()?;
    module.add_class::<PyHarmonicSimilarityOptions>()?;
    module.add_class::<PyHarmonicSimilarity>()?;
    module.add_class::<PyGroupVariance>()?;
    module.add_class::<PyConvergenceBlock>()?;
    module.add_function(wrap_pyfunction!(analyse_rmsd_to_reference, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_pairwise_fitted_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_pairwise_torus_distance, module)?)?;
    module.add_function(wrap_pyfunction!(
        analyse_generalized_procrustes_mean,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(analyse_agglomerative_clustering, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_dbscan_clustering, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_medoid, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_kmeans, module)?)?;
    module.add_function(wrap_pyfunction!(
        analyse_harmonic_ensemble_similarity,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        analyse_cluster_population_similarity,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(analyse_group_coordinate_variance, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_block_convergence, module)?)?;
    Ok(())
}
