//! Owned Python projection of native MRC density maps and statistics.

use crate::crystallography::PyUnitCell;
use numpy::ndarray::ArrayView1;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyclass(name = "MapBoundary", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyMapBoundary {
    Missing,
    Periodic,
}

#[pyclass(name = "MapStatistics", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMapStatistics {
    #[pyo3(get)]
    pub(crate) count: usize,
    #[pyo3(get)]
    pub(crate) minimum: f32,
    #[pyo3(get)]
    pub(crate) maximum: f32,
    #[pyo3(get)]
    pub(crate) mean: f64,
    #[pyo3(get)]
    pub(crate) sigma: f64,
}

#[pyclass(name = "MapHistogram", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMapHistogram {
    #[pyo3(get)]
    pub(crate) minimum: f32,
    #[pyo3(get)]
    pub(crate) maximum: f32,
    #[pyo3(get)]
    pub(crate) counts: Vec<usize>,
}

#[pyclass(name = "DensityMap", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDensityMap(pub(crate) molframe::xtal::DensityMap);

#[pymethods]
impl PyDensityMap {
    #[new]
    #[pyo3(signature = (dimensions, starts, sampling, cell, origin, space_group, labels, extended_header, values))]
    fn new(
        dimensions: [usize; 3],
        starts: [i32; 3],
        sampling: [usize; 3],
        cell: &PyUnitCell,
        origin: [f64; 3],
        space_group: i32,
        labels: Vec<String>,
        extended_header: Vec<u8>,
        values: Vec<f32>,
    ) -> Self {
        Self(molframe::xtal::DensityMap {
            dimensions,
            starts,
            sampling,
            cell: cell.cell,
            origin,
            space_group,
            labels: labels.into_iter().map(Into::into).collect(),
            extended_header,
            values,
        })
    }

    #[staticmethod]
    fn from_mrc_bytes(data: &[u8]) -> PyResult<Self> {
        molframe::xtal::DensityMap::from_mrc_bytes(data)
            .map(Self)
            .map_err(value_error)
    }

    #[getter]
    fn dimensions(&self) -> [usize; 3] {
        self.0.dimensions
    }

    #[getter]
    fn starts(&self) -> [i32; 3] {
        self.0.starts
    }

    #[getter]
    fn sampling(&self) -> [usize; 3] {
        self.0.sampling
    }

    #[getter]
    fn cell(&self) -> PyResult<PyUnitCell> {
        PyUnitCell::from_native(self.0.cell)
    }

    #[getter]
    fn origin(&self) -> [f64; 3] {
        self.0.origin
    }

    #[getter]
    fn space_group(&self) -> i32 {
        self.0.space_group
    }

    #[getter]
    fn labels(&self) -> Vec<String> {
        self.0.labels.iter().map(ToString::to_string).collect()
    }

    #[getter]
    fn extended_header(&self) -> Vec<u8> {
        self.0.extended_header.clone()
    }

    #[getter]
    fn values<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        let (length, pointer) = {
            let density_map = slf.borrow();
            (density_map.0.values.len(), density_map.0.values.as_ptr())
        };
        // SAFETY: the frozen density map owns a contiguous f32 vector and is
        // retained as the NumPy base below for the complete view lifetime.
        let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
        // SAFETY: the map cannot mutate its vector after construction and the
        // returned array is made read-only before it reaches Python.
        let values = unsafe { PyArray1::borrow_from_array(&view, slf.clone().into_any()) };
        let _readonly = values.readwrite().make_nonwriteable();
        values
    }

    fn to_mrc_bytes(&self, py: Python<'_>) -> PyResult<Py<PyBytes>> {
        self.0
            .to_mrc_bytes()
            .map(|data| PyBytes::new(py, &data).unbind())
            .map_err(value_error)
    }

    fn statistics(&self) -> PyResult<PyMapStatistics> {
        self.0.statistics().map(Into::into).map_err(value_error)
    }

    fn masked_statistics(&self, mask: Vec<bool>) -> PyResult<PyMapStatistics> {
        self.0
            .masked_statistics(&mask)
            .map(Into::into)
            .map_err(value_error)
    }

    fn histogram(&self, bins: usize, minimum: f32, maximum: f32) -> PyResult<PyMapHistogram> {
        self.0
            .histogram(bins, minimum, maximum)
            .map(Into::into)
            .map_err(value_error)
    }
}

impl From<molframe::xtal::MapStatistics> for PyMapStatistics {
    fn from(value: molframe::xtal::MapStatistics) -> Self {
        Self {
            count: value.count,
            minimum: value.minimum,
            maximum: value.maximum,
            mean: value.mean,
            sigma: value.sigma,
        }
    }
}

impl From<molframe::xtal::MapHistogram> for PyMapHistogram {
    fn from(value: molframe::xtal::MapHistogram) -> Self {
        Self {
            minimum: value.minimum,
            maximum: value.maximum,
            counts: value.counts,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMapBoundary>()?;
    module.add_class::<PyMapStatistics>()?;
    module.add_class::<PyMapHistogram>()?;
    module.add_class::<PyDensityMap>()?;
    Ok(())
}
