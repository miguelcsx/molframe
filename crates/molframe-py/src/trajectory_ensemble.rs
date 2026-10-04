//! Ensemble analyses of a trajectory held as an array: distances, principal components,
//! clusters, embeddings, displacement and convergence, each with its provenance.

// pyo3 hands the functions below owned lists extracted from Python sequences.
#![allow(clippy::needless_pass_by_value)]

mod results;
mod statistics;

use self::results::{
    PyClustering, PyConvergence, PyDiffusionMap, PyGroupVariance, PyHarmonicSimilarity, PyKMeans,
    PyMsd, PyPathSimilarity, PyPca, matrix,
};
use crate::analysis_result::PyAnalysis;
use crate::policy::{PyAnalysisPolicy, policy_of};
use molframe::trajectory::{
    CartesianFit, EnsembleDistanceMatrix, FrameView, KMeansOptions, Linkage,
    analyse_agglomerative_clustering, analyse_cartesian_pca_view, analyse_dbscan_clustering,
    analyse_diffusion_map, analyse_generalized_procrustes_mean_view, analyse_kmeans_view,
    analyse_medoid, analyse_pairwise_fitted_rmsd_view,
};
use numpy::{IntoPyArray, PyArray2, PyArray3, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;
use std::sync::Arc;

/// Default ceiling for the dense matrices an ensemble analysis allocates: 1 GiB.
pub(super) const MEMORY_LIMIT: usize = 1 << 30;

pub(super) fn frames<'a>(array: &'a numpy::PyReadonlyArray3<'_, f32>) -> PyResult<FrameView<'a>> {
    let shape = array.shape();
    let flat = array.as_slice().map_err(|_| {
        crate::error::value("positions must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    if shape[2] != 3 {
        return Err(crate::error::value(
            "positions must have shape (frames, atoms, 3)",
        ));
    }
    FrameView::new(flat.as_chunks::<3>().0, shape[0], shape[1]).map_err(crate::error::kernel)
}

fn distances(array: &Bound<'_, PyArray2<f64>>) -> PyResult<EnsembleDistanceMatrix> {
    let readonly = array.readonly();
    let view = readonly.as_array();
    if view.nrows() != view.ncols() {
        return Err(crate::error::value("distances must be a square matrix"));
    }
    Ok(EnsembleDistanceMatrix {
        size: view.nrows(),
        values: view.iter().copied().collect(),
    })
}

pub(super) fn wrap<T>(
    py: Python<'_>,
    analysis: &molframe::Analysis<T>,
    value: impl FnOnce(Python<'_>, &T) -> PyResult<Py<PyAny>>,
) -> PyResult<PyAnalysis> {
    let converted = match analysis.value() {
        Some(inner) => Some(value(py, inner)?),
        None => None,
    };
    Ok(PyAnalysis::new(analysis, converted))
}

pub(super) fn rows_of(array: &Bound<'_, PyArray2<f64>>) -> Vec<Vec<f64>> {
    let readonly = array.readonly();
    readonly
        .as_array()
        .rows()
        .into_iter()
        .map(|row| row.to_vec())
        .collect()
}

/// Pairwise RMSD between every two frames after the best rigid fit of each pair.
///
/// `positions` has shape `(frames, atoms, 3)`. The matrix is `(frames, frames)`;
/// `memory_limit` is the ceiling in bytes on what it may allocate.
#[pyfunction]
#[pyo3(signature = (positions, *, memory_limit=MEMORY_LIMIT, policy=None))]
fn pairwise_rmsd(
    py: Python<'_>,
    positions: &Bound<'_, PyArray3<f32>>,
    memory_limit: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = positions.readonly();
    let view = frames(&array)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_pairwise_fitted_rmsd_view(view, memory_limit, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(matrix(py, value.size, value.size, value.values.to_vec())?
            .into_any()
            .unbind())
    })
}

/// Principal components of the Cartesian coordinates of every frame.
///
/// `fit` chooses how frames are aligned first: `"none"` keeps them as given, `"mean"` fits to
/// an iterative Procrustes mean (`tolerance`, `max_iterations`), and a `reference` array of
/// shape `(atoms, 3)` fits to that structure.
#[pyfunction]
#[pyo3(signature = (
    positions,
    *,
    components,
    fit="none",
    reference=None,
    tolerance=1e-6,
    max_iterations=100,
    memory_limit=MEMORY_LIMIT,
    policy=None,
))]
#[allow(clippy::too_many_arguments)]
fn pca(
    py: Python<'_>,
    positions: &Bound<'_, PyArray3<f32>>,
    components: usize,
    fit: &str,
    reference: Option<&Bound<'_, PyArray2<f32>>>,
    tolerance: f64,
    max_iterations: usize,
    memory_limit: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = positions.readonly();
    let view = frames(&array)?;
    let reference = reference.map(numpy::PyArrayMethods::readonly);
    let reference = match &reference {
        Some(reference) => Some(crate::bindings::coordinates(reference)?),
        None => None,
    };
    let fit = match (fit, reference) {
        ("none", None) => CartesianFit::None,
        ("mean", None) => CartesianFit::IterativeMean {
            tolerance,
            max_iterations,
        },
        (_, Some(reference)) => CartesianFit::Reference(reference),
        (other, None) => {
            return Err(crate::error::value(format!(
                "{other:?} is not a fit; the fits are \"none\" and \"mean\", or pass a reference"
            )));
        }
    };
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_cartesian_pca_view(view, fit, components, memory_limit, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(Py::new(
            py,
            PyPca {
                inner: Arc::new(value.clone()),
            },
        )?
        .into_any())
    })
}

/// The mean structure of the frames, from generalized Procrustes alignment.
#[pyfunction]
#[pyo3(signature = (positions, *, tolerance=1e-6, max_iterations=100, policy=None))]
fn mean_structure(
    py: Python<'_>,
    positions: &Bound<'_, PyArray3<f32>>,
    tolerance: f64,
    max_iterations: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = positions.readonly();
    let view = frames(&array)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| {
            analyse_generalized_procrustes_mean_view(view, tolerance, max_iterations, &policy)
        })
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        let flat: Vec<f32> = value.iter().flatten().copied().collect();
        let array = flat.into_pyarray(py).reshape((value.len(), 3))?;
        array.readwrite().make_nonwriteable();
        Ok(array.into_any().unbind())
    })
}

/// A diffusion map of a square distance matrix: affinities `exp(-d^2 / epsilon)`.
#[pyfunction]
#[pyo3(signature = (distances, *, epsilon, dimensions, time=1, policy=None))]
fn diffusion_map(
    py: Python<'_>,
    distances: &Bound<'_, PyArray2<f64>>,
    epsilon: f64,
    dimensions: usize,
    time: u32,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let matrix = self::distances(distances)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_diffusion_map(&matrix, "supplied", epsilon, time, dimensions, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(Py::new(
            py,
            PyDiffusionMap {
                inner: Arc::new(value.clone()),
            },
        )?
        .into_any())
    })
}

/// K-means on a `(observations, features)` matrix, started from the observations listed in
/// `initial_centres` (one per cluster, so the start is a decision the caller made).
#[pyfunction]
#[pyo3(signature = (
    observations,
    *,
    initial_centres,
    maximum_iterations=300,
    tolerance_squared=1e-12,
    policy=None,
))]
fn kmeans(
    py: Python<'_>,
    observations: &Bound<'_, PyArray2<f64>>,
    initial_centres: Vec<usize>,
    maximum_iterations: usize,
    tolerance_squared: f64,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let readonly = observations.readonly();
    let shape = readonly.shape();
    let (rows, columns) = (shape[0], shape[1]);
    let values = readonly.as_slice().map_err(|_| {
        crate::error::value("observations must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    let options = KMeansOptions {
        initial_centres: &initial_centres,
        maximum_iterations,
        convergence_tolerance_squared: tolerance_squared,
    };
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_kmeans_view(values, rows, columns, options, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(Py::new(
            py,
            PyKMeans {
                inner: Arc::new(value.clone()),
            },
        )?
        .into_any())
    })
}

fn clustering(py: Python<'_>, value: &molframe::trajectory::Clustering) -> PyResult<Py<PyAny>> {
    Ok(Py::new(
        py,
        PyClustering {
            inner: Arc::new(value.clone()),
        },
    )?
    .into_any())
}

/// Agglomerative clustering of a square distance matrix into `clusters` clusters.
///
/// `linkage` is `"single"`, `"complete"` or `"average"`.
#[pyfunction]
#[pyo3(signature = (distances, *, clusters, linkage="average", policy=None))]
fn agglomerative_clustering(
    py: Python<'_>,
    distances: &Bound<'_, PyArray2<f64>>,
    clusters: usize,
    linkage: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let matrix = self::distances(distances)?;
    let linkage: Linkage = linkage.parse().map_err(crate::error::kernel)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_agglomerative_clustering(&matrix, clusters, linkage, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, clustering)
}

/// DBSCAN on a square distance matrix: points within `epsilon` of at least `minimum_points`
/// others seed a cluster; the rest are noise (label `-1`).
#[pyfunction]
#[pyo3(signature = (distances, *, epsilon, minimum_points, policy=None))]
fn dbscan(
    py: Python<'_>,
    distances: &Bound<'_, PyArray2<f64>>,
    epsilon: f64,
    minimum_points: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let matrix = self::distances(distances)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_dbscan_clustering(&matrix, epsilon, minimum_points, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, clustering)
}

/// The member with the least total distance to the other members.
#[pyfunction]
#[pyo3(signature = (distances, members, *, policy=None))]
fn medoid(
    py: Python<'_>,
    distances: &Bound<'_, PyArray2<f64>>,
    members: Vec<usize>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let matrix = self::distances(distances)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_medoid(&matrix, &members, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(value.into_pyobject(py)?.into_any().unbind())
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPca>()?;
    module.add_class::<PyDiffusionMap>()?;
    module.add_class::<PyKMeans>()?;
    module.add_class::<PyClustering>()?;
    module.add_class::<PyMsd>()?;
    module.add_class::<PyConvergence>()?;
    module.add_class::<PyGroupVariance>()?;
    module.add_class::<PyHarmonicSimilarity>()?;
    module.add_class::<PyPathSimilarity>()?;
    module.add_function(wrap_pyfunction!(pairwise_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(pca, module)?)?;
    module.add_function(wrap_pyfunction!(mean_structure, module)?)?;
    module.add_function(wrap_pyfunction!(diffusion_map, module)?)?;
    module.add_function(wrap_pyfunction!(kmeans, module)?)?;
    module.add_function(wrap_pyfunction!(agglomerative_clustering, module)?)?;
    module.add_function(wrap_pyfunction!(dbscan, module)?)?;
    module.add_function(wrap_pyfunction!(medoid, module)?)?;
    statistics::register(module)
}
