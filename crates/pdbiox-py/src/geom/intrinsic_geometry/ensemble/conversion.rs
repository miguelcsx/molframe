//! Native-to-Python ensemble result conversion.

use super::{
    PyClustering, PyConvergenceBlock, PyEnsembleDistanceMatrix, PyFrameAlignment, PyGroupVariance,
    PyHarmonicSimilarity, PyHarmonicSimilarityOptions, PyKMeans, PyKMeansOptions, PyLinkage,
    PyRemainderPolicy,
};
use pyo3::IntoPyObjectExt;
use pyo3::prelude::{Py, PyAny, PyResult, Python};
use pyo3::types::PyList;

impl PyEnsembleDistanceMatrix {
    pub(crate) fn native(&self) -> pdbiox::traj::EnsembleDistanceMatrix {
        pdbiox::traj::EnsembleDistanceMatrix {
            size: self.size,
            values: self.values.clone().into_boxed_slice(),
        }
    }
}
impl PyKMeansOptions {
    pub(crate) fn native(&self) -> pdbiox::traj::KMeansOptions<'_> {
        pdbiox::traj::KMeansOptions {
            initial_centres: &self.initial_centres,
            maximum_iterations: self.maximum_iterations,
            convergence_tolerance_squared: self.convergence_tolerance_squared,
        }
    }
}
impl PyHarmonicSimilarityOptions {
    pub(crate) fn native(&self) -> pdbiox::traj::HarmonicSimilarityOptions {
        self.0
    }
}
impl From<PyFrameAlignment> for pdbiox::traj::FrameAlignment {
    fn from(v: PyFrameAlignment) -> Self {
        match v {
            PyFrameAlignment::Unaligned => Self::None,
            PyFrameAlignment::Rigid => Self::Rigid,
        }
    }
}
impl From<PyLinkage> for pdbiox::traj::Linkage {
    fn from(v: PyLinkage) -> Self {
        match v {
            PyLinkage::Single => Self::Single,
            PyLinkage::Complete => Self::Complete,
            PyLinkage::Average => Self::Average,
        }
    }
}
impl From<PyRemainderPolicy> for pdbiox::traj::RemainderPolicy {
    fn from(v: PyRemainderPolicy) -> Self {
        match v {
            PyRemainderPolicy::Include => Self::Include,
            PyRemainderPolicy::Reject => Self::Reject,
        }
    }
}
impl From<pdbiox::traj::Clustering> for PyClustering {
    fn from(v: pdbiox::traj::Clustering) -> Self {
        Self {
            labels: v.labels,
            members: v.members,
            medoids: v.medoids,
        }
    }
}
impl From<pdbiox::traj::KMeans> for PyKMeans {
    fn from(v: pdbiox::traj::KMeans) -> Self {
        Self {
            labels: v.labels,
            centres: v.centres,
            inertia: v.inertia,
            iterations: v.iterations,
        }
    }
}
impl From<pdbiox::traj::EnsembleDistanceMatrix> for PyEnsembleDistanceMatrix {
    fn from(value: pdbiox::traj::EnsembleDistanceMatrix) -> Self {
        Self {
            size: value.size,
            values: value.values.into_vec(),
        }
    }
}
impl From<pdbiox::traj::HarmonicSimilarity> for PyHarmonicSimilarity {
    fn from(v: pdbiox::traj::HarmonicSimilarity) -> Self {
        Self {
            mean_term: v.mean_term,
            covariance_term: v.covariance_term,
            distance: v.distance,
            similarity: v.similarity,
        }
    }
}
impl From<pdbiox::traj::GroupVariance> for PyGroupVariance {
    fn from(v: pdbiox::traj::GroupVariance) -> Self {
        Self {
            atoms: v.atoms,
            variance_by_axis: v.variance_by_axis,
            rms_fluctuation: v.rms_fluctuation,
        }
    }
}
impl From<pdbiox::traj::ConvergenceBlock> for PyConvergenceBlock {
    fn from(v: pdbiox::traj::ConvergenceBlock) -> Self {
        Self {
            start: v.start,
            end: v.end,
            block_mean: v.block_mean,
            cumulative_mean: v.cumulative_mean,
        }
    }
}
pub(super) fn distance_value(
    py: Python<'_>,
    v: pdbiox::traj::EnsembleDistanceMatrix,
) -> PyResult<Py<PyAny>> {
    Py::new(
        py,
        PyEnsembleDistanceMatrix {
            size: v.size,
            values: v.values.into_vec(),
        },
    )
    .map(Py::into_any)
}
pub(super) fn clustering_value(
    py: Python<'_>,
    value: pdbiox::traj::Clustering,
) -> PyResult<Py<PyAny>> {
    Py::new(py, PyClustering::from(value)).map(Py::into_any)
}
pub(super) fn kmeans_value(py: Python<'_>, value: pdbiox::traj::KMeans) -> PyResult<Py<PyAny>> {
    Py::new(py, PyKMeans::from(value)).map(Py::into_any)
}
pub(super) fn harmonic_value(
    py: Python<'_>,
    value: pdbiox::traj::HarmonicSimilarity,
) -> PyResult<Py<PyAny>> {
    Py::new(py, PyHarmonicSimilarity::from(value)).map(Py::into_any)
}
pub(super) fn float_list_value(py: Python<'_>, values: Vec<f64>) -> PyResult<Py<PyAny>> {
    Ok(PyList::new(py, values)?.unbind().into_any())
}
pub(super) fn usize_value(py: Python<'_>, value: usize) -> PyResult<Py<PyAny>> {
    value.into_py_any(py)
}
pub(super) fn f64_value(py: Python<'_>, value: f64) -> PyResult<Py<PyAny>> {
    value.into_py_any(py)
}
pub(super) fn coordinate_list_value(py: Python<'_>, values: Vec<[f32; 3]>) -> PyResult<Py<PyAny>> {
    Ok(PyList::new(py, values)?.unbind().into_any())
}
pub(super) fn group_list_value(
    py: Python<'_>,
    values: Vec<pdbiox::traj::GroupVariance>,
) -> PyResult<Py<PyAny>> {
    let values = values
        .into_iter()
        .map(|v| Py::new(py, PyGroupVariance::from(v)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, values)?.unbind().into_any())
}
pub(super) fn block_list_value(
    py: Python<'_>,
    values: Vec<pdbiox::traj::ConvergenceBlock>,
) -> PyResult<Py<PyAny>> {
    let values = values
        .into_iter()
        .map(|v| Py::new(py, PyConvergenceBlock::from(v)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, values)?.unbind().into_any())
}
