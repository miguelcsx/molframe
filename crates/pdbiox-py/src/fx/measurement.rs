//! Native constraint-measurement values and execution.

use super::alignment::PyAlignedMotif;
use super::errors::measurement_error;
use super::specification::PyMotif;
use crate::geometry::PyEigenOptions;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "MeasurementValue", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMeasurementValue(pub(crate) pdbiox::fx::MeasurementValue);

#[pymethods]
impl PyMeasurementValue {
    #[staticmethod]
    fn scalar(value: f64) -> PyResult<Self> {
        if !value.is_finite() {
            return Err(PyValueError::new_err("scalar measurement must be finite"));
        }
        Ok(Self(pdbiox::fx::MeasurementValue::Scalar(value)))
    }

    #[staticmethod]
    fn chirality(sign: i8) -> PyResult<Self> {
        if !matches!(sign, -1 | 1) {
            return Err(PyValueError::new_err("chirality sign must be -1 or 1"));
        }
        Ok(Self(pdbiox::fx::MeasurementValue::Chirality(sign)))
    }

    #[staticmethod]
    fn count(value: usize) -> Self {
        Self(pdbiox::fx::MeasurementValue::Count(value))
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.0 {
            pdbiox::fx::MeasurementValue::Scalar(_) => "scalar",
            pdbiox::fx::MeasurementValue::Chirality(_) => "chirality",
            pdbiox::fx::MeasurementValue::Count(_) => "count",
            _ => "unknown",
        }
    }

    #[getter]
    fn numeric(&self) -> f64 {
        self.0.numeric()
    }

    #[getter]
    fn scalar_value(&self) -> Option<f64> {
        match self.0 {
            pdbiox::fx::MeasurementValue::Scalar(value) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn chirality_sign(&self) -> Option<i8> {
        match self.0 {
            pdbiox::fx::MeasurementValue::Chirality(value) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn count_value(&self) -> Option<usize> {
        match self.0 {
            pdbiox::fx::MeasurementValue::Count(value) => Some(value),
            _ => None,
        }
    }
}

#[pyclass(name = "IndeterminateReason", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyIndeterminateReason {
    MissingCoordinate,
    DegenerateGeometry,
}

impl From<pdbiox::fx::IndeterminateReason> for PyIndeterminateReason {
    fn from(value: pdbiox::fx::IndeterminateReason) -> Self {
        match value {
            pdbiox::fx::IndeterminateReason::MissingCoordinate => Self::MissingCoordinate,
            pdbiox::fx::IndeterminateReason::DegenerateGeometry => Self::DegenerateGeometry,
        }
    }
}

#[pyclass(name = "ConstraintMeasurement", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyConstraintMeasurement {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    value: Option<PyMeasurementValue>,
    #[pyo3(get)]
    deviation: Option<f64>,
    #[pyo3(get)]
    satisfied: Option<bool>,
    #[pyo3(get)]
    atoms: Vec<u32>,
    #[pyo3(get)]
    indeterminate: Option<PyIndeterminateReason>,
}

impl From<pdbiox::fx::ConstraintMeasurement> for PyConstraintMeasurement {
    fn from(value: pdbiox::fx::ConstraintMeasurement) -> Self {
        Self {
            name: value.name.to_string(),
            value: value.value.map(PyMeasurementValue),
            deviation: value.deviation,
            satisfied: value.satisfied,
            atoms: value
                .atoms
                .into_iter()
                .map(pdbiox::AtomIndex::get)
                .collect(),
            indeterminate: value.indeterminate.map(Into::into),
        }
    }
}

#[pyclass(name = "MeasurementSet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMeasurementSet(pub(crate) pdbiox::fx::MeasurementSet);

#[pymethods]
impl PyMeasurementSet {
    #[getter]
    fn constraints(&self) -> Vec<PyConstraintMeasurement> {
        self.0.constraints.iter().cloned().map(Into::into).collect()
    }

    #[getter]
    fn metrics(&self) -> BTreeMap<String, f64> {
        self.0
            .metrics()
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect()
    }
}

#[pyclass(name = "MeasurementOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMeasurementOptions(pub(crate) pdbiox::fx::MeasurementOptions);

#[pymethods]
impl PyMeasurementOptions {
    #[new]
    #[pyo3(signature = (maximum_alternatives, *, plane_fit=None))]
    fn new(maximum_alternatives: usize, plane_fit: Option<&PyEigenOptions>) -> PyResult<Self> {
        if maximum_alternatives == 0 {
            return Err(PyValueError::new_err(
                "maximum_alternatives must be greater than zero",
            ));
        }
        Ok(Self(pdbiox::fx::MeasurementOptions {
            maximum_alternatives,
            plane_fit: plane_fit.map_or_else(pdbiox::EigenOptions::standard, |value| value.inner),
        }))
    }

    #[getter]
    const fn maximum_alternatives(&self) -> usize {
        self.0.maximum_alternatives
    }

    #[getter]
    fn plane_fit(&self) -> PyEigenOptions {
        PyEigenOptions {
            inner: self.0.plane_fit,
        }
    }
}

#[pyfunction]
pub(crate) fn measure_constraints(
    py: Python<'_>,
    structure: &PyStructure,
    motif: &PyMotif,
    aligned: &PyAlignedMotif,
    options: &PyMeasurementOptions,
) -> PyResult<PyMeasurementSet> {
    let structure = structure.structure().clone();
    let motif = motif.0.clone();
    let aligned = aligned.0.clone();
    let options = options.0;
    py.detach(move || {
        pdbiox::fx::measure_constraints(&structure, &motif, &aligned, options)
            .map(PyMeasurementSet)
            .map_err(measurement_error)
    })
}
