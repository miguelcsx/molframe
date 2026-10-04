//! Displacement, variance, convergence and similarity of ensembles.

use super::results::{
    PyConvergence, PyGroupVariance, PyHarmonicSimilarity, PyMsd, PyPathSimilarity,
};
use super::{MEMORY_LIMIT, frames, rows_of, wrap};
use crate::analysis_result::PyAnalysis;
use crate::policy::{PyAnalysisPolicy, policy_of};
use molframe::trajectory::{
    HarmonicSimilarityOptions, PathFrameMetric, RemainderPolicy, analyse_block_convergence,
    analyse_cluster_population_similarity, analyse_group_coordinate_variance_view,
    analyse_harmonic_ensemble_similarity, analyse_mean_squared_displacement_view,
    path_similarity_view,
};
use numpy::{PyArray2, PyArray3, PyArrayMethods};
use pyo3::prelude::*;
use std::sync::Arc;

/// Window-averaged mean squared displacement for lags `0..=maximum_lag`.
///
/// Coordinates are used as given: unwrap a periodic trajectory first. `atoms` selects atoms
/// by index; the default is every atom.
#[pyfunction]
#[pyo3(signature = (positions, *, maximum_lag, atoms=None, policy=None))]
fn msd(
    py: Python<'_>,
    positions: &Bound<'_, PyArray3<f32>>,
    maximum_lag: usize,
    atoms: Option<Vec<usize>>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = positions.readonly();
    let view = frames(&array)?;
    let atoms: Vec<usize> = atoms.into_iter().flatten().collect();
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_mean_squared_displacement_view(view, &atoms, maximum_lag, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(Py::new(
            py,
            PyMsd {
                rows: Arc::new(value.clone()),
            },
        )?
        .into_any())
    })
}

/// Coordinate variance of each atom group over frames that are already aligned.
///
/// Frames are not fitted here: variance only means something relative to the alignment the
/// caller chose.
#[pyfunction]
#[pyo3(signature = (positions, groups, *, policy=None))]
fn group_variance(
    py: Python<'_>,
    positions: &Bound<'_, PyArray3<f32>>,
    groups: Vec<Vec<usize>>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = positions.readonly();
    let view = frames(&array)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_group_coordinate_variance_view(view, &groups, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        let found: Vec<Py<PyGroupVariance>> = value
            .iter()
            .map(|inner| {
                Py::new(
                    py,
                    PyGroupVariance {
                        inner: inner.clone(),
                    },
                )
            })
            .collect::<PyResult<_>>()?;
        Ok(found.into_pyobject(py)?.into_any().unbind())
    })
}

/// Block means and cumulative means of a series, to see whether it has converged.
///
/// `remainder` is `"include"` (a final partial block, with its exact bounds) or `"reject"`.
#[pyfunction]
#[pyo3(signature = (values, *, block_size, remainder="include", policy=None))]
fn block_convergence(
    py: Python<'_>,
    values: Vec<f64>,
    block_size: usize,
    remainder: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let remainder: RemainderPolicy = remainder.parse().map_err(crate::error::kernel)?;
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_block_convergence(&values, block_size, remainder, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(Py::new(
            py,
            PyConvergence {
                blocks: Arc::new(value.clone()),
            },
        )?
        .into_any())
    })
}

/// Harmonic ensemble similarity: Gaussian models fitted to two sets of observations.
///
/// `regularization` is added to every covariance diagonal element.
#[pyfunction]
#[pyo3(signature = (first, second, *, regularization, memory_limit=MEMORY_LIMIT, policy=None))]
fn harmonic_similarity(
    py: Python<'_>,
    first: &Bound<'_, PyArray2<f64>>,
    second: &Bound<'_, PyArray2<f64>>,
    regularization: f64,
    memory_limit: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let (first, second) = (rows_of(first), rows_of(second));
    let options = HarmonicSimilarityOptions {
        covariance_regularization: regularization,
        memory_limit_bytes: memory_limit,
    };
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_harmonic_ensemble_similarity(&first, &second, options, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(Py::new(py, PyHarmonicSimilarity { inner: *value })?.into_any())
    })
}

/// Jensen-Shannon similarity of the cluster populations of two ensembles.
///
/// Cluster identities must already be aligned. One means identical, zero disjoint.
#[pyfunction]
#[pyo3(signature = (first, second, *, clusters, policy=None))]
fn population_similarity(
    py: Python<'_>,
    first: Vec<usize>,
    second: Vec<usize>,
    clusters: usize,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let policy = policy_of(policy);
    let analysis = py
        .detach(|| analyse_cluster_population_similarity(&first, &second, clusters, &policy))
        .map_err(crate::error::kernel)?;
    wrap(py, &analysis, |py, value| {
        Ok(value.into_pyobject(py)?.into_any().unbind())
    })
}

/// Hausdorff and discrete Fréchet distances between two paths of frames.
///
/// `metric` is `"fitted_rmsd"` (each pair of frames rigidly fitted first) or
/// `"cartesian_rmsd"`.
#[pyfunction]
#[pyo3(signature = (first, second, *, metric="fitted_rmsd", memory_limit=MEMORY_LIMIT))]
fn path_similarity(
    py: Python<'_>,
    first: &Bound<'_, PyArray3<f32>>,
    second: &Bound<'_, PyArray3<f32>>,
    metric: &str,
    memory_limit: usize,
) -> PyResult<PyPathSimilarity> {
    let metric: PathFrameMetric = metric.parse().map_err(crate::error::kernel)?;
    let (first, second) = (first.readonly(), second.readonly());
    let (first, second) = (frames(&first)?, frames(&second)?);
    let inner = py
        .detach(|| path_similarity_view(first, second, metric, memory_limit))
        .map_err(crate::error::kernel)?;
    Ok(PyPathSimilarity { inner })
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(msd, module)?)?;
    module.add_function(wrap_pyfunction!(group_variance, module)?)?;
    module.add_function(wrap_pyfunction!(block_convergence, module)?)?;
    module.add_function(wrap_pyfunction!(harmonic_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(population_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(path_similarity, module)?)
}
