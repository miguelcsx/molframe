//! Coarse Python calls into native dynamics kernels.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::{PyAnalysisPolicy, PySelection};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "WaterDynamicsOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyWaterDynamicsOptions {
    #[pyo3(get)]
    maximum_lag: usize,
}

#[pymethods]
impl PyWaterDynamicsOptions {
    #[new]
    fn new(maximum_lag: usize) -> Self {
        Self { maximum_lag }
    }
}

impl PyWaterDynamicsOptions {
    const fn native(self) -> molframe::analysis::WaterDynamicsOptions {
        molframe::analysis::WaterDynamicsOptions {
            maximum_lag: self.maximum_lag,
        }
    }
}

#[pyclass(name = "WaterLag", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyWaterLag {
    #[pyo3(get)]
    lag: usize,
    #[pyo3(get)]
    observations: u64,
    #[pyo3(get)]
    surviving: u64,
    #[pyo3(get)]
    resident: u64,
    #[pyo3(get)]
    survival_probability: Option<f64>,
    #[pyo3(get)]
    residence_probability: Option<f64>,
}

#[pyclass(name = "DielectricOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDielectricOptions(pub(crate) molframe::analysis::DielectricOptions);

#[pymethods]
impl PyDielectricOptions {
    #[new]
    fn new(volume: f64, temperature: f64, fluctuation_prefactor: f64) -> Self {
        Self(molframe::analysis::DielectricOptions {
            volume,
            temperature,
            fluctuation_prefactor,
        })
    }

    #[getter]
    fn volume(&self) -> f64 {
        self.0.volume
    }

    #[getter]
    fn temperature(&self) -> f64 {
        self.0.temperature
    }

    #[getter]
    fn fluctuation_prefactor(&self) -> f64 {
        self.0.fluctuation_prefactor
    }
}

#[pyclass(name = "DielectricResult", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDielectricResult {
    #[pyo3(get)]
    mean_dipole: [f64; 3],
    #[pyo3(get)]
    fluctuation: f64,
    #[pyo3(get)]
    relative_permittivity: f64,
}

#[pyfunction]
pub(crate) fn water_dynamics(
    py: Python<'_>,
    occupancy: Vec<PySelection>,
    options: PyWaterDynamicsOptions,
) -> PyResult<Vec<PyWaterLag>> {
    let occupancy: Vec<_> = occupancy.into_iter().map(|value| value.inner).collect();
    py.detach(move || molframe::analysis::water_dynamics(&occupancy, options.native()))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_water_dynamics(
    py: Python<'_>,
    occupancy: Vec<PySelection>,
    options: PyWaterDynamicsOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let occupancy: Vec<_> = occupancy.into_iter().map(|value| value.inner).collect();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            molframe::analysis::governed_water_dynamics(&occupancy, options.native(), &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, water_lag_list)
}

#[pyfunction]
pub(crate) fn dielectric_from_dipoles(
    py: Python<'_>,
    dipoles: Vec<[f64; 3]>,
    options: PyDielectricOptions,
) -> PyResult<PyDielectricResult> {
    py.detach(move || molframe::analysis::dielectric_from_dipoles(&dipoles, options.0))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_dielectric_from_dipoles(
    py: Python<'_>,
    dipoles: Vec<[f64; 3]>,
    options: PyDielectricOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            molframe::analysis::governed_dielectric_from_dipoles(&dipoles, options.0, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Ok(Py::new(py, PyDielectricResult::from(value))?.into_any())
    })
}

impl From<molframe::analysis::WaterLag> for PyWaterLag {
    fn from(value: molframe::analysis::WaterLag) -> Self {
        Self {
            lag: value.lag,
            observations: value.observations,
            surviving: value.surviving,
            resident: value.resident,
            survival_probability: value.survival_probability,
            residence_probability: value.residence_probability,
        }
    }
}

impl From<molframe::analysis::DielectricResult> for PyDielectricResult {
    fn from(value: molframe::analysis::DielectricResult) -> Self {
        Self {
            mean_dipole: value.mean_dipole,
            fluctuation: value.fluctuation,
            relative_permittivity: value.relative_permittivity,
        }
    }
}

fn water_lag_list(
    py: Python<'_>,
    values: Vec<molframe::analysis::WaterLag>,
) -> PyResult<Py<PyAny>> {
    let values = values
        .into_iter()
        .map(|value| Py::new(py, PyWaterLag::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, values)?.unbind().into_any())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
