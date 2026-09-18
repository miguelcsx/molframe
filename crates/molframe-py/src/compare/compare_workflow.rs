//! Explicit comparison workflows delegated to the facade's native kernels.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::errors::{compare_error, governed_compare_error};
use crate::geometry::{PyRigid, borrowed_coordinates};
use crate::query::PyAnalysisPolicy;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

#[pyclass(name = "PointMatch", frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyPointMatch {
    #[pyo3(get)]
    reference: usize,
    #[pyo3(get)]
    model: usize,
}

#[pymethods]
impl PyPointMatch {
    #[new]
    const fn new(reference: usize, model: usize) -> Self {
        Self { reference, model }
    }
}

#[pyclass(name = "PointMapping", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPointMapping {
    inner: molframe::compare::PointMapping,
    reference_len: usize,
    model_len: usize,
}

#[pymethods]
impl PyPointMapping {
    #[new]
    fn new(matches: Vec<PyPointMatch>, reference_len: usize, model_len: usize) -> PyResult<Self> {
        let native = matches.into_iter().map(native_match);
        let inner = molframe::compare::PointMapping::new(native, reference_len, model_len)
            .map_err(|error| compare_error(&error))?;
        Ok(Self {
            inner,
            reference_len,
            model_len,
        })
    }

    #[getter]
    fn matches(&self) -> Vec<PyPointMatch> {
        self.inner
            .matches()
            .iter()
            .copied()
            .map(Into::into)
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.matches().len()
    }
}

impl PyPointMapping {
    pub(crate) fn inner(&self) -> &molframe::compare::PointMapping {
        &self.inner
    }

    fn validate_lengths(&self, reference_len: usize, model_len: usize) -> PyResult<()> {
        if self.reference_len != reference_len || self.model_len != model_len {
            return Err(PyValueError::new_err(format!(
                "mapping was created for ({}, {}) points, received ({}, {})",
                self.reference_len, self.model_len, reference_len, model_len
            )));
        }
        Ok(())
    }
}

#[pyclass(name = "ComparisonAlignment", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyComparisonAlignment(molframe::compare::ComparisonAlignment);

#[pymethods]
impl PyComparisonAlignment {
    #[staticmethod]
    fn not_required() -> Self {
        Self(molframe::compare::ComparisonAlignment::NotRequired)
    }

    #[staticmethod]
    fn rigid(transform: PyRigid) -> Self {
        Self(molframe::compare::ComparisonAlignment::Rigid(transform.0))
    }

    #[getter]
    fn is_required(&self) -> bool {
        matches!(self.0, molframe::compare::ComparisonAlignment::Rigid(_))
    }

    #[getter]
    fn transform(&self) -> Option<PyRigid> {
        match self.0 {
            molframe::compare::ComparisonAlignment::NotRequired => None,
            molframe::compare::ComparisonAlignment::Rigid(transform) => Some(PyRigid(transform)),
        }
    }
}

impl PyComparisonAlignment {
    fn inner(self) -> molframe::compare::ComparisonAlignment {
        self.0
    }
}

#[pyclass(name = "DistanceMeasurement", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDistanceMeasurement {
    #[pyo3(get)]
    distances: Vec<f64>,
    #[pyo3(get)]
    rmsd: f64,
}

#[pymethods]
impl PyDistanceMeasurement {
    #[new]
    fn new(distances: Vec<f64>, rmsd: f64) -> Self {
        Self { distances, rmsd }
    }
}

impl From<molframe::compare::DistanceMeasurement> for PyDistanceMeasurement {
    fn from(value: molframe::compare::DistanceMeasurement) -> Self {
        Self {
            distances: value.distances.into_vec(),
            rmsd: value.rmsd,
        }
    }
}

#[pyclass(name = "ComparisonVerdict", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyComparisonVerdict {
    #[pyo3(get)]
    maximum_rmsd: f64,
    #[pyo3(get)]
    passed: bool,
}

#[pymethods]
impl PyComparisonVerdict {
    #[new]
    const fn new(maximum_rmsd: f64, passed: bool) -> Self {
        Self {
            maximum_rmsd,
            passed,
        }
    }
}

impl From<molframe::compare::ComparisonVerdict> for PyComparisonVerdict {
    fn from(value: molframe::compare::ComparisonVerdict) -> Self {
        Self {
            maximum_rmsd: value.maximum_rmsd,
            passed: value.passed,
        }
    }
}

#[pyclass(name = "InterfaceRmsd", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyInterfaceRmsd {
    #[pyo3(get)]
    value: f64,
}

#[pymethods]
impl PyInterfaceRmsd {
    #[new]
    const fn new(value: f64) -> Self {
        Self { value }
    }
}

#[pyclass(name = "PocketRmsd", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPocketRmsd {
    #[pyo3(get)]
    value: f64,
}

#[pymethods]
impl PyPocketRmsd {
    #[new]
    const fn new(value: f64) -> Self {
        Self { value }
    }
}

#[pyfunction]
pub(crate) fn align_mapping(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
) -> PyResult<PyComparisonAlignment> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    py.detach(|| molframe::compare::align_mapping(reference, model, mapping.inner()))
        .map(Into::into)
        .map_err(|error| compare_error(&error))
}

#[pyfunction]
pub(crate) fn measure_mapping(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
    alignment: &PyComparisonAlignment,
) -> PyResult<PyDistanceMeasurement> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    let alignment = alignment.inner();
    Ok(py
        .detach(|| molframe::compare::measure_mapping(reference, model, mapping.inner(), alignment))
        .into())
}

#[pyfunction]
pub(crate) fn decide_rmsd(
    py: Python<'_>,
    measurement: &PyDistanceMeasurement,
    maximum_rmsd: f64,
) -> PyComparisonVerdict {
    py.detach(move || -> PyComparisonVerdict {
        let native = molframe::compare::DistanceMeasurement {
            distances: measurement.distances.clone().into_boxed_slice(),
            rmsd: measurement.rmsd,
        };
        molframe::compare::decide_rmsd(&native, maximum_rmsd).into()
    })
}

#[pyfunction]
pub(crate) fn interface_rmsd(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
    alignment: &PyComparisonAlignment,
) -> PyResult<PyInterfaceRmsd> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    let alignment = alignment.inner();
    Ok(py
        .detach(|| molframe::compare::interface_rmsd(reference, model, mapping.inner(), alignment))
        .into())
}

#[pyfunction]
pub(crate) fn pocket_rmsd(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
    alignment: &PyComparisonAlignment,
) -> PyResult<PyPocketRmsd> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    let alignment = alignment.inner();
    Ok(py
        .detach(|| molframe::compare::pocket_rmsd(reference, model, mapping.inner(), alignment))
        .into())
}

#[pyfunction]
pub(crate) fn governed_comparison_workflow(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
    maximum_rmsd: f64,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            molframe::compare::governed_comparison_workflow(
                reference,
                model,
                mapping.inner(),
                maximum_rmsd,
                &policy,
            )
        })
        .map_err(|error| governed_compare_error(&error))?;
    analysis_with_value(py, result, |py, (measurement, verdict)| {
        let measurement = Py::new(py, PyDistanceMeasurement::from(measurement))?;
        let verdict = Py::new(py, PyComparisonVerdict::from(verdict))?;
        Ok(
            PyTuple::new(py, [measurement.into_any(), verdict.into_any()])?
                .unbind()
                .into_any(),
        )
    })
}

#[pyfunction]
pub(crate) fn governed_interface_rmsd(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            molframe::compare::governed_interface_rmsd(reference, model, mapping.inner(), &policy)
        })
        .map_err(|error| governed_compare_error(&error))?;
    analysis_with_value(py, result, |py, value| {
        Py::new(py, PyInterfaceRmsd::from(value)).map(Py::into_any)
    })
}

#[pyfunction]
pub(crate) fn governed_pocket_rmsd(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    model: PyReadonlyArray2<'_, f32>,
    mapping: &PyPointMapping,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = borrowed_coordinates(&reference)?;
    let model = borrowed_coordinates(&model)?;
    mapping.validate_lengths(reference.len(), model.len())?;
    let policy = policy.inner.clone();
    let result = py
        .detach(|| {
            molframe::compare::governed_pocket_rmsd(reference, model, mapping.inner(), &policy)
        })
        .map_err(|error| governed_compare_error(&error))?;
    analysis_with_value(py, result, |py, value| {
        Py::new(py, PyPocketRmsd::from(value)).map(Py::into_any)
    })
}

fn native_match(value: PyPointMatch) -> molframe::compare::PointMatch {
    molframe::compare::PointMatch {
        reference: value.reference,
        model: value.model,
    }
}

impl From<molframe::compare::PointMatch> for PyPointMatch {
    fn from(value: molframe::compare::PointMatch) -> Self {
        Self {
            reference: value.reference,
            model: value.model,
        }
    }
}

impl From<molframe::compare::ComparisonAlignment> for PyComparisonAlignment {
    fn from(value: molframe::compare::ComparisonAlignment) -> Self {
        Self(value)
    }
}

impl From<molframe::compare::InterfaceRmsd> for PyInterfaceRmsd {
    fn from(value: molframe::compare::InterfaceRmsd) -> Self {
        Self { value: value.0 }
    }
}

impl From<molframe::compare::PocketRmsd> for PyPocketRmsd {
    fn from(value: molframe::compare::PocketRmsd) -> Self {
        Self { value: value.0 }
    }
}
