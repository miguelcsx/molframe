//! Native trajectory result projection without Python data-plane loops.

use super::model::PyTrajectoryOperation;
use crate::contract::{PyAnalysis, analysis_with_value};
use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Columnar `NumPy` representation of a mean-squared-displacement series.
#[pyclass(name = "MeanSquaredDisplacementSeries", frozen)]
pub(crate) struct PyMeanSquaredDisplacementSeries {
    lags: Py<PyArray1<usize>>,
    observations: Py<PyArray1<u64>>,
    values: Py<PyArray1<f64>>,
}

#[pymethods]
impl PyMeanSquaredDisplacementSeries {
    #[getter]
    fn lags(&self, py: Python<'_>) -> Py<PyArray1<usize>> {
        self.lags.clone_ref(py)
    }

    #[getter]
    fn observations(&self, py: Python<'_>) -> Py<PyArray1<u64>> {
        self.observations.clone_ref(py)
    }

    #[getter]
    fn values(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.values.clone_ref(py)
    }
}

impl PyMeanSquaredDisplacementSeries {
    fn from_native(py: Python<'_>, values: Vec<pdbiox::traj::MeanSquaredDisplacement>) -> Self {
        let mut lags = Vec::with_capacity(values.len());
        let mut observations = Vec::with_capacity(values.len());
        let mut series = Vec::with_capacity(values.len());
        for value in values {
            lags.push(value.lag);
            observations.push(value.observations);
            series.push(value.value);
        }
        Self {
            lags: Array1::from_vec(lags).into_pyarray(py).unbind(),
            observations: Array1::from_vec(observations).into_pyarray(py).unbind(),
            values: Array1::from_vec(series).into_pyarray(py).unbind(),
        }
    }
}

pub(crate) fn execute_operation(
    py: Python<'_>,
    operation: &PyTrajectoryOperation,
) -> PyResult<PyAnalysis> {
    let plan = super::super::PyPlan {
        operations: vec![(
            "result".to_owned(),
            super::super::Operation::Trajectory(operation.clone_ref(py)),
        )],
    };
    let native = super::super::execute_native(&plan, py, None, None)?;
    let Some(entry) = native.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native trajectory plan returned no result",
        ));
    };
    match entry.value {
        pdbiox::PlanValue::Trajectory(value) => analysis_to_python(py, *value),
        _ => Err(PyValueError::new_err(
            "native trajectory plan returned an incompatible result",
        )),
    }
}

pub(crate) fn value_to_python(
    py: Python<'_>,
    value: pdbiox::TrajectoryValue,
) -> PyResult<Py<PyAny>> {
    Py::new(py, analysis_to_python(py, value)?).map(Py::into_any)
}

fn analysis_to_python(py: Python<'_>, value: pdbiox::TrajectoryValue) -> PyResult<PyAnalysis> {
    match value {
        pdbiox::TrajectoryValue::RmsdToReference(analysis) => {
            analysis_with_value(py, analysis, |py, values| Ok(float_vector(py, values)))
        }
        pdbiox::TrajectoryValue::MeanSquaredDisplacement(analysis) => {
            analysis_with_value(py, analysis, msd_series)
        }
        pdbiox::TrajectoryValue::PairwiseFittedRmsd(analysis) => {
            analysis_with_value(py, analysis, dense_matrix)
        }
        pdbiox::TrajectoryValue::GeneralizedProcrustesMean(analysis) => {
            analysis_with_value(py, analysis, coordinate_matrix)
        }
    }
}

fn float_vector(py: Python<'_>, values: Vec<f64>) -> Py<PyAny> {
    Array1::from_vec(values)
        .into_pyarray(py)
        .unbind()
        .into_any()
}

fn msd_series(
    py: Python<'_>,
    values: Vec<pdbiox::traj::MeanSquaredDisplacement>,
) -> PyResult<Py<PyAny>> {
    Py::new(py, PyMeanSquaredDisplacementSeries::from_native(py, values)).map(Py::into_any)
}

fn dense_matrix(
    py: Python<'_>,
    values: pdbiox::traj::EnsembleDistanceMatrix,
) -> PyResult<Py<PyAny>> {
    Array2::from_shape_vec((values.size, values.size), values.values.into_vec())
        .map(|values| values.into_pyarray(py).unbind().into_any())
        .map_err(|_| PyValueError::new_err("native pairwise matrix has an invalid shape"))
}

fn coordinate_matrix(py: Python<'_>, values: Vec<[f32; 3]>) -> PyResult<Py<PyAny>> {
    let capacity = values.len().checked_mul(3).ok_or_else(|| {
        PyValueError::new_err("Procrustes result exceeds addressable NumPy shape")
    })?;
    let mut flat = Vec::with_capacity(capacity);
    for coordinate in values {
        flat.extend(coordinate);
    }
    Array2::from_shape_vec((capacity / 3, 3), flat)
        .map(|values| values.into_pyarray(py).unbind().into_any())
        .map_err(|_| PyValueError::new_err("native Procrustes result has an invalid shape"))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMeanSquaredDisplacementSeries>()
}
