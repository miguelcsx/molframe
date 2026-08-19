//! Coarse, vectorised Python projections of governed ensemble kernels.

use super::{borrowed_frames, value_error};
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::PyAnalysisPolicy;
use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3, PyUntypedArrayMethods};
use pdbiox::{PeriodicAngle, TorusMetric};
use pyo3::prelude::*;

#[pyclass(name = "FrameAlignment", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyFrameAlignment {
    Unaligned,
    Rigid,
}

#[pyclass(name = "Linkage", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyLinkage {
    Single,
    Complete,
    Average,
}

#[pyclass(name = "RemainderPolicy", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyRemainderPolicy {
    Include,
    Reject,
}

#[pyclass(name = "EnsembleDistanceMatrix", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEnsembleDistanceMatrix {
    #[pyo3(get)]
    size: usize,
    #[pyo3(get)]
    values: Vec<f64>,
}

#[pymethods]
impl PyEnsembleDistanceMatrix {
    #[new]
    fn new(values: PyReadonlyArray2<'_, f64>) -> PyResult<Self> {
        let shape = values.shape();
        if shape[0] != shape[1] {
            return Err(value_error("distance matrix must be square"));
        }
        let result = Self {
            size: shape[0],
            values: values.as_array().iter().copied().collect(),
        };
        drop(values);
        Ok(result)
    }
}

#[pyclass(name = "Clustering", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyClustering {
    #[pyo3(get)]
    labels: Vec<Option<usize>>,
    #[pyo3(get)]
    members: Vec<Vec<usize>>,
    #[pyo3(get)]
    medoids: Vec<usize>,
}

#[pyclass(name = "KMeansOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyKMeansOptions {
    initial_centres: Vec<usize>,
    maximum_iterations: usize,
    convergence_tolerance_squared: f64,
}

#[pymethods]
impl PyKMeansOptions {
    #[new]
    fn new(
        initial_centres: Vec<usize>,
        maximum_iterations: usize,
        convergence_tolerance_squared: f64,
    ) -> Self {
        Self {
            initial_centres,
            maximum_iterations,
            convergence_tolerance_squared,
        }
    }
}

#[pyclass(name = "KMeans", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyKMeans {
    #[pyo3(get)]
    labels: Vec<usize>,
    #[pyo3(get)]
    centres: Vec<Vec<f64>>,
    #[pyo3(get)]
    inertia: f64,
    #[pyo3(get)]
    iterations: usize,
}

#[pyclass(name = "HarmonicSimilarityOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHarmonicSimilarityOptions(pdbiox::traj::HarmonicSimilarityOptions);

#[pymethods]
impl PyHarmonicSimilarityOptions {
    #[new]
    fn new(covariance_regularization: f64, memory_limit_bytes: usize) -> Self {
        Self(pdbiox::traj::HarmonicSimilarityOptions {
            covariance_regularization,
            memory_limit_bytes,
        })
    }
}

#[pyclass(name = "HarmonicSimilarity", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHarmonicSimilarity {
    #[pyo3(get)]
    mean_term: f64,
    #[pyo3(get)]
    covariance_term: f64,
    #[pyo3(get)]
    distance: f64,
    #[pyo3(get)]
    similarity: f64,
}

#[pyclass(name = "GroupVariance", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGroupVariance {
    #[pyo3(get)]
    atoms: Vec<usize>,
    #[pyo3(get)]
    variance_by_axis: [f64; 3],
    #[pyo3(get)]
    rms_fluctuation: f64,
}

#[pyclass(name = "ConvergenceBlock", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyConvergenceBlock {
    #[pyo3(get)]
    start: usize,
    #[pyo3(get)]
    end: usize,
    #[pyo3(get)]
    block_mean: f64,
    #[pyo3(get)]
    cumulative_mean: f64,
}

#[pyfunction]
pub(crate) fn analyse_rmsd_to_reference(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    reference: usize,
    alignment: PyFrameAlignment,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let frames = borrowed_frames(&frames)?;
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            pdbiox::traj::analyse_rmsd_to_reference_view(
                frames,
                reference,
                alignment.into(),
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, float_list_value)
}

#[pyfunction]
pub(crate) fn analyse_pairwise_fitted_rmsd(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    memory_limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let frames = borrowed_frames(&frames)?;
    let policy = policy.inner.clone();
    let result = py
        .detach(|| pdbiox::traj::analyse_pairwise_fitted_rmsd_view(frames, memory_limit, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, result, distance_value)
}

#[pyfunction]
pub(crate) fn analyse_pairwise_torus_distance(
    py: Python<'_>,
    points: PyReadonlyArray2<'_, f64>,
    weights: PyReadonlyArray1<'_, f64>,
    metric_name: &str,
    memory_limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let point_values = points
        .as_array()
        .outer_iter()
        .map(|row| {
            row.iter()
                .copied()
                .map(PeriodicAngle::from_radians)
                .collect()
        })
        .collect::<Result<Vec<Vec<_>>, _>>()
        .map_err(value_error)?;
    let metric = TorusMetric::new(weights.as_array().iter().copied().collect::<Vec<_>>())
        .map_err(value_error)?;
    drop((points, weights));
    let metric_name = metric_name.to_owned();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            pdbiox::traj::analyse_pairwise_torus_distance(
                &point_values,
                &metric,
                &metric_name,
                memory_limit,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, distance_value)
}

#[pyfunction]
pub(crate) fn analyse_generalized_procrustes_mean(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    tolerance: f64,
    maximum_iterations: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let frames = borrowed_frames(&frames)?;
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            pdbiox::traj::analyse_generalized_procrustes_mean_view(
                frames,
                tolerance,
                maximum_iterations,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, coordinate_list_value)
}

#[pyfunction]
pub(crate) fn analyse_agglomerative_clustering(
    py: Python<'_>,
    distances: &PyEnsembleDistanceMatrix,
    cluster_count: usize,
    linkage: PyLinkage,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let distances = distances.native();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            pdbiox::traj::analyse_agglomerative_clustering(
                &distances,
                cluster_count,
                linkage.into(),
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, clustering_value)
}

#[pyfunction]
pub(crate) fn analyse_dbscan_clustering(
    py: Python<'_>,
    distances: &PyEnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let distances = distances.native();
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            pdbiox::traj::analyse_dbscan_clustering(&distances, epsilon, minimum_points, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, clustering_value)
}

#[pyfunction]
pub(crate) fn analyse_medoid(
    py: Python<'_>,
    distances: &PyEnsembleDistanceMatrix,
    members: Vec<usize>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let distances = distances.native();
    let policy = policy.inner.clone();
    let result = py
        .detach(move || pdbiox::traj::analyse_medoid(&distances, &members, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, result, usize_value)
}

#[pyfunction]
pub(crate) fn analyse_kmeans(
    py: Python<'_>,
    observations: PyReadonlyArray2<'_, f64>,
    options: &PyKMeansOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let shape = observations.shape();
    let rows = shape[0];
    let columns = shape[1];
    let observations = observations.as_slice().map_err(|_| {
        value_error("trajectory arrays must be C-contiguous; pass an explicit contiguous copy")
    })?;
    let options = options.clone();
    let policy = policy.inner.clone();
    let result = py
        .detach(move || {
            pdbiox::traj::analyse_kmeans_view(
                observations,
                rows,
                columns,
                options.native(),
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, kmeans_value)
}

#[pyfunction]
pub(crate) fn analyse_harmonic_ensemble_similarity(
    py: Python<'_>,
    first: Vec<Vec<f64>>,
    second: Vec<Vec<f64>>,
    options: PyHarmonicSimilarityOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let policy = policy.inner.clone();
    let result = py
        .detach(move || {
            pdbiox::traj::analyse_harmonic_ensemble_similarity(&first, &second, options.0, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, harmonic_value)
}

#[pyfunction]
pub(crate) fn analyse_cluster_population_similarity(
    py: Python<'_>,
    first: Vec<usize>,
    second: Vec<usize>,
    cluster_count: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let policy = policy.inner.clone();
    let result = py
        .detach(move || {
            pdbiox::traj::analyse_cluster_population_similarity(
                &first,
                &second,
                cluster_count,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, f64_value)
}

#[pyfunction]
pub(crate) fn analyse_group_coordinate_variance(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    groups: Vec<Vec<usize>>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let frames = borrowed_frames(&frames)?;
    let policy = policy.inner.clone();
    let result = py
        .detach(move || {
            pdbiox::traj::analyse_group_coordinate_variance_view(frames, &groups, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, group_list_value)
}

#[pyfunction]
pub(crate) fn analyse_block_convergence(
    py: Python<'_>,
    values: Vec<f64>,
    block_size: usize,
    remainder: PyRemainderPolicy,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let policy = policy.inner.clone();
    let result = py
        .detach(move || {
            pdbiox::traj::analyse_block_convergence(&values, block_size, remainder.into(), &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, result, block_list_value)
}
#[path = "ensemble/conversion.rs"]
mod conversion;
use conversion::{
    block_list_value, clustering_value, coordinate_list_value, distance_value, f64_value,
    float_list_value, group_list_value, harmonic_value, kmeans_value, usize_value,
};
