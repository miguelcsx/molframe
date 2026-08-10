//! Typed Python projection of calibration-free helical geometry.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::PyAnalysisPolicy;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "BaseFrame", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBaseFrame(pub(crate) pdbiox::analysis::BaseFrame);

#[pymethods]
impl PyBaseFrame {
    #[new]
    fn new(origin: [f64; 3], x: [f64; 3], y: [f64; 3], z: [f64; 3]) -> Self {
        Self(pdbiox::analysis::BaseFrame { origin, x, y, z })
    }

    #[getter]
    fn origin(&self) -> [f64; 3] {
        self.0.origin
    }
    #[getter]
    fn x(&self) -> [f64; 3] {
        self.0.x
    }
    #[getter]
    fn y(&self) -> [f64; 3] {
        self.0.y
    }
    #[getter]
    fn z(&self) -> [f64; 3] {
        self.0.z
    }
}

#[pyclass(name = "HelicalOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHelicalOptions {
    #[pyo3(get)]
    frame_tolerance: f64,
}

#[pymethods]
impl PyHelicalOptions {
    #[new]
    fn new(frame_tolerance: f64) -> Self {
        Self { frame_tolerance }
    }
}

impl PyHelicalOptions {
    const fn native(self) -> pdbiox::analysis::HelicalOptions {
        pdbiox::analysis::HelicalOptions {
            frame_tolerance: self.frame_tolerance,
        }
    }
}

#[pyclass(name = "HelicalParameters", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHelicalParameters {
    #[pyo3(get)]
    x_displacement: f64,
    #[pyo3(get)]
    y_displacement: f64,
    #[pyo3(get)]
    z_displacement: f64,
    #[pyo3(get)]
    x_rotation_degrees: f64,
    #[pyo3(get)]
    y_rotation_degrees: f64,
    #[pyo3(get)]
    z_rotation_degrees: f64,
}

#[pyfunction]
pub(crate) fn helical_parameters(
    py: Python<'_>,
    first: PyBaseFrame,
    second: PyBaseFrame,
    options: PyHelicalOptions,
) -> PyResult<PyHelicalParameters> {
    py.detach(move || pdbiox::analysis::helical_parameters(first.0, second.0, options.native()))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn helical_steps(
    py: Python<'_>,
    frames: Vec<PyBaseFrame>,
    options: PyHelicalOptions,
) -> PyResult<Vec<PyHelicalParameters>> {
    let frames: Vec<_> = frames.into_iter().map(|value| value.0).collect();
    py.detach(move || pdbiox::analysis::helical_steps(&frames, options.native()))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_helical_parameters(
    py: Python<'_>,
    first: PyBaseFrame,
    second: PyBaseFrame,
    options: PyHelicalOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            pdbiox::analysis::governed_helical_parameters(
                first.0,
                second.0,
                options.native(),
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, parameter_value)
}

#[pyfunction]
pub(crate) fn analyse_helical_steps(
    py: Python<'_>,
    frames: Vec<PyBaseFrame>,
    options: PyHelicalOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let frames: Vec<_> = frames.into_iter().map(|value| value.0).collect();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            pdbiox::analysis::governed_helical_steps(&frames, options.native(), &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, parameter_list)
}

impl From<pdbiox::analysis::HelicalParameters> for PyHelicalParameters {
    fn from(value: pdbiox::analysis::HelicalParameters) -> Self {
        Self {
            x_displacement: value.x_displacement,
            y_displacement: value.y_displacement,
            z_displacement: value.z_displacement,
            x_rotation_degrees: value.x_rotation_degrees,
            y_rotation_degrees: value.y_rotation_degrees,
            z_rotation_degrees: value.z_rotation_degrees,
        }
    }
}

fn parameter_value(
    py: Python<'_>,
    value: pdbiox::analysis::HelicalParameters,
) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyHelicalParameters::from(value))?.into_any())
}

fn parameter_list(
    py: Python<'_>,
    values: Vec<pdbiox::analysis::HelicalParameters>,
) -> PyResult<Py<PyAny>> {
    let values = values
        .into_iter()
        .map(|value| Py::new(py, PyHelicalParameters::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, values)?.unbind().into_any())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
