//! Python projection of `molframe_core::contract::Analysis`.

use super::{PyDiagnostic, PyProvenance};
use molframe::{Analysis, Assumption, AssumptionSource, Coverage, ImpactEstimate, Status};
use pyo3::prelude::*;

#[pyclass(name = "Status", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyStatus {
    Complete,
    Partial,
    Ambiguous,
    Indeterminate,
}

impl From<Status> for PyStatus {
    fn from(value: Status) -> Self {
        match value {
            Status::Complete => Self::Complete,
            Status::Partial => Self::Partial,
            Status::Ambiguous => Self::Ambiguous,
            _ => Self::Indeterminate,
        }
    }
}

#[pyclass(name = "Coverage", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCoverage {
    #[pyo3(get)]
    intended: u32,
    #[pyo3(get)]
    used: u32,
    #[pyo3(get)]
    missing: u32,
    #[pyo3(get)]
    ambiguous: u32,
}

#[pymethods]
impl PyCoverage {
    #[getter]
    fn fraction(&self) -> f64 {
        Coverage::from(*self).fraction()
    }

    #[getter]
    fn is_complete(&self) -> bool {
        Coverage::from(*self).is_complete()
    }
}

impl From<Coverage> for PyCoverage {
    fn from(value: Coverage) -> Self {
        Self {
            intended: value.intended,
            used: value.used,
            missing: value.missing,
            ambiguous: value.ambiguous,
        }
    }
}

impl From<PyCoverage> for Coverage {
    fn from(value: PyCoverage) -> Self {
        Self {
            intended: value.intended,
            used: value.used,
            missing: value.missing,
            ambiguous: value.ambiguous,
        }
    }
}

#[pyclass(name = "ImpactEstimate", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyImpactEstimate {
    NoImpact,
    Low,
    Moderate,
    High,
    Unknown,
}

impl From<ImpactEstimate> for PyImpactEstimate {
    fn from(value: ImpactEstimate) -> Self {
        match value {
            ImpactEstimate::None => Self::NoImpact,
            ImpactEstimate::Low => Self::Low,
            ImpactEstimate::Moderate => Self::Moderate,
            ImpactEstimate::High => Self::High,
            _ => Self::Unknown,
        }
    }
}

#[pyclass(name = "AssumptionSource", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAssumptionSource {
    Explicit,
    ProfileDefault,
    Inferred,
}

impl From<AssumptionSource> for PyAssumptionSource {
    fn from(value: AssumptionSource) -> Self {
        match value {
            AssumptionSource::Explicit => Self::Explicit,
            AssumptionSource::ProfileDefault => Self::ProfileDefault,
            _ => Self::Inferred,
        }
    }
}

#[pyclass(name = "Assumption", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAssumption {
    #[pyo3(get)]
    field: String,
    #[pyo3(get)]
    value: String,
    #[pyo3(get)]
    source: PyAssumptionSource,
    #[pyo3(get)]
    impact: PyImpactEstimate,
    #[pyo3(get)]
    is_silent_and_material: bool,
}

impl From<Assumption> for PyAssumption {
    fn from(value: Assumption) -> Self {
        let is_silent_and_material = value.is_silent_and_material();
        Self {
            field: value.field.name().to_owned(),
            value: value.value.into(),
            source: value.source.into(),
            impact: value.impact.into(),
            is_silent_and_material,
        }
    }
}

#[pyclass(name = "Analysis", frozen)]
pub(crate) struct PyAnalysis {
    value: Py<PyAny>,
    #[pyo3(get)]
    status: PyStatus,
    #[pyo3(get)]
    coverage: PyCoverage,
    #[pyo3(get)]
    warnings: Vec<PyDiagnostic>,
    #[pyo3(get)]
    assumptions: Vec<PyAssumption>,
    #[pyo3(get)]
    provenance: PyProvenance,
}

#[pymethods]
impl PyAnalysis {
    #[getter]
    fn value(&self, py: Python<'_>) -> Py<PyAny> {
        self.value.clone_ref(py)
    }

    #[getter]
    fn is_usable(&self) -> bool {
        self.status != PyStatus::Indeterminate
    }

    #[getter]
    fn silent_assumptions(&self) -> Vec<PyAssumption> {
        self.assumptions
            .iter()
            .filter(|item| item.is_silent_and_material)
            .cloned()
            .collect()
    }
}

pub(crate) fn analysis_with_value<T>(
    py: Python<'_>,
    analysis: Analysis<T>,
    convert: impl FnOnce(Python<'_>, T) -> PyResult<Py<PyAny>>,
) -> PyResult<PyAnalysis> {
    Ok(PyAnalysis {
        value: convert(py, analysis.value)?,
        status: analysis.status.into(),
        coverage: analysis.coverage.into(),
        warnings: analysis
            .warnings
            .into_iter()
            .map(PyDiagnostic::from)
            .collect(),
        assumptions: analysis
            .assumptions
            .into_iter()
            .map(PyAssumption::from)
            .collect(),
        provenance: analysis.provenance.into(),
    })
}
