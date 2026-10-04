//! The result types of the ensemble analyses: each holds its arrays and hands them out
//! read-only.

use molframe::trajectory::{ConvergenceBlock, GroupVariance, HarmonicSimilarity, PathSimilarity};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use std::sync::Arc;

pub(super) fn matrix(
    py: Python<'_>,
    rows: usize,
    columns: usize,
    flat: Vec<f64>,
) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let array = flat.into_pyarray(py).reshape((rows, columns))?;
    array.readwrite().make_nonwriteable();
    Ok(array)
}

/// A `(rows, columns)` array from equal-length rows.
fn table<'py, R: AsRef<[f64]>>(py: Python<'py>, rows: &[R]) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let columns = rows.first().map_or(0, |row| row.as_ref().len());
    let flat = rows
        .iter()
        .flat_map(|row| row.as_ref().iter().copied())
        .collect();
    matrix(py, rows.len(), columns, flat)
}

fn vector<T: numpy::Element>(py: Python<'_>, values: Vec<T>) -> Bound<'_, PyArray1<T>> {
    let array = values.into_pyarray(py);
    array.readwrite().make_nonwriteable();
    array
}

/// Principal components of a set of observations.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Pca",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyPca {
    pub(super) inner: Arc<molframe::trajectory::PcaResult>,
}

#[pymethods]
impl PyPca {
    /// The mean of every feature.
    #[getter]
    fn mean<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vector(py, self.inner.mean.to_vec())
    }

    /// Eigenvalues, largest first.
    #[getter]
    fn eigenvalues<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vector(py, self.inner.eigenvalues.to_vec())
    }

    /// One axis per component, `(components, features)`, each with a canonical sign.
    #[getter]
    fn components<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        table(py, &self.inner.components)
    }

    /// Each observation's coordinates on the components, `(observations, components)`.
    #[getter]
    fn projections<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        table(py, &self.inner.projections)
    }
}

/// A diffusion-map embedding of an affinity geometry.
#[derive(Clone, Debug)]
#[pyclass(
    name = "DiffusionMap",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyDiffusionMap {
    pub(super) inner: Arc<molframe::trajectory::DiffusionMap>,
}

#[pymethods]
impl PyDiffusionMap {
    /// Non-constant eigenvalues, largest first.
    #[getter]
    fn eigenvalues<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vector(py, self.inner.eigenvalues.to_vec())
    }

    /// Observation coordinates after the diffusion time, `(observations, dimensions)`.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        table(py, &self.inner.coordinates)
    }

    /// The connected component of the affinity graph each observation lies in.
    #[getter]
    fn graph_components<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        vector(
            py,
            self.inner
                .graph_components
                .iter()
                .map(|&c| c as u64)
                .collect(),
        )
    }
}

/// A k-means partition.
#[derive(Clone, Debug)]
#[pyclass(
    name = "KMeans",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyKMeans {
    pub(super) inner: Arc<molframe::trajectory::KMeans>,
}

#[pymethods]
impl PyKMeans {
    /// The cluster of each observation.
    #[getter]
    fn labels<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        vector(
            py,
            self.inner
                .labels
                .iter()
                .map(|&label| label as u64)
                .collect(),
        )
    }

    /// The mean of each cluster, `(clusters, features)`.
    #[getter]
    fn centres<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        table(py, &self.inner.centres)
    }

    /// Sum of squared distances to the assigned centres.
    #[getter]
    fn inertia(&self) -> f64 {
        self.inner.inertia
    }

    /// Iterations needed to converge.
    #[getter]
    fn iterations(&self) -> usize {
        self.inner.iterations
    }
}

/// A partition into clusters, each with a medoid representative.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Clustering",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyClustering {
    pub(super) inner: Arc<molframe::trajectory::Clustering>,
}

#[pymethods]
impl PyClustering {
    /// The cluster of each observation; `-1` marks DBSCAN noise.
    #[getter]
    fn labels<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i64>> {
        let labels = self
            .inner
            .labels
            .iter()
            .map(|label| match label {
                Some(cluster) => match i64::try_from(*cluster) {
                    Ok(cluster) => cluster,
                    Err(_) => i64::MAX,
                },
                None => -1,
            })
            .collect();
        vector(py, labels)
    }

    /// The sorted observation indices of each cluster.
    #[getter]
    fn members(&self) -> Vec<Vec<usize>> {
        self.inner.members.clone()
    }

    /// The member of each cluster with the least total distance to the others.
    #[getter]
    fn medoids(&self) -> Vec<usize> {
        self.inner.medoids.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.members.len()
    }
}

/// Mean squared displacement by frame lag.
#[derive(Clone, Debug)]
#[pyclass(
    name = "MeanSquaredDisplacement",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyMsd {
    pub(super) rows: Arc<Vec<molframe::trajectory::MeanSquaredDisplacement>>,
}

#[pymethods]
impl PyMsd {
    /// The frame lag of each row.
    #[getter]
    fn lag<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        vector(py, self.rows.iter().map(|row| row.lag as u64).collect())
    }

    /// How many atom-origin observations each row averages.
    #[getter]
    fn observations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        vector(py, self.rows.iter().map(|row| row.observations).collect())
    }

    /// Mean squared Cartesian displacement at each lag.
    #[getter]
    fn value<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vector(py, self.rows.iter().map(|row| row.value).collect())
    }

    fn __len__(&self) -> usize {
        self.rows.len()
    }
}

/// Mean per-block and cumulative values of a series.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Convergence",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyConvergence {
    pub(super) blocks: Arc<Vec<ConvergenceBlock>>,
}

#[pymethods]
impl PyConvergence {
    /// First observation of each block (inclusive).
    #[getter]
    fn start<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        vector(
            py,
            self.blocks.iter().map(|block| block.start as u64).collect(),
        )
    }

    /// Last observation of each block (exclusive).
    #[getter]
    fn end<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        vector(
            py,
            self.blocks.iter().map(|block| block.end as u64).collect(),
        )
    }

    /// The mean within each block.
    #[getter]
    fn block_mean<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vector(
            py,
            self.blocks.iter().map(|block| block.block_mean).collect(),
        )
    }

    /// The mean from the first observation through each block.
    #[getter]
    fn cumulative_mean<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vector(
            py,
            self.blocks
                .iter()
                .map(|block| block.cumulative_mean)
                .collect(),
        )
    }

    fn __len__(&self) -> usize {
        self.blocks.len()
    }
}

/// The coordinate variance of one atom group.
#[derive(Clone, Debug)]
#[pyclass(
    name = "GroupVariance",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyGroupVariance {
    pub(super) inner: GroupVariance,
}

#[pymethods]
impl PyGroupVariance {
    /// The atoms of the group, in the order given.
    #[getter]
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms.clone()
    }

    /// Mean population variance along x, y and z over the group's atoms.
    #[getter]
    fn variance_by_axis(&self) -> [f64; 3] {
        self.inner.variance_by_axis
    }

    /// Square root of the sum of the three axis variances.
    #[getter]
    fn rms_fluctuation(&self) -> f64 {
        self.inner.rms_fluctuation
    }
}

/// A distance between two Gaussian models of an ensemble.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "HarmonicSimilarity",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyHarmonicSimilarity {
    pub(super) inner: HarmonicSimilarity,
}

#[pymethods]
impl PyHarmonicSimilarity {
    /// The mean-separation contribution to the distance.
    #[getter]
    const fn mean_term(&self) -> f64 {
        self.inner.mean_term
    }

    /// The covariance-volume contribution to the distance.
    #[getter]
    const fn covariance_term(&self) -> f64 {
        self.inner.covariance_term
    }

    /// The sum of the two terms (a Bhattacharyya distance).
    #[getter]
    const fn distance(&self) -> f64 {
        self.inner.distance
    }

    /// `exp(-distance)`: one for identical models.
    #[getter]
    const fn similarity(&self) -> f64 {
        self.inner.similarity
    }
}

/// How far apart two paths are.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "PathSimilarity",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(super) struct PyPathSimilarity {
    pub(super) inner: PathSimilarity,
}

#[pymethods]
impl PyPathSimilarity {
    #[getter]
    const fn hausdorff_distance(&self) -> f64 {
        self.inner.hausdorff_distance
    }

    #[getter]
    const fn discrete_frechet_distance(&self) -> f64 {
        self.inner.discrete_frechet_distance
    }
}
