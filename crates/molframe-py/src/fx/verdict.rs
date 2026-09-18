//! Declarative verdict profiles and their decomposed outcomes.

use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "Comparison", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyComparison(pub(crate) molframe::fx::Comparison);

#[pymethods]
impl PyComparison {
    #[staticmethod]
    fn less_than(limit: f64) -> PyResult<Self> {
        finite_limit(limit).map(|()| Self(molframe::fx::Comparison::LessThan(limit)))
    }

    #[staticmethod]
    fn at_most(limit: f64) -> PyResult<Self> {
        finite_limit(limit).map(|()| Self(molframe::fx::Comparison::AtMost(limit)))
    }

    #[staticmethod]
    fn greater_than(limit: f64) -> PyResult<Self> {
        finite_limit(limit).map(|()| Self(molframe::fx::Comparison::GreaterThan(limit)))
    }

    #[staticmethod]
    fn at_least(limit: f64) -> PyResult<Self> {
        finite_limit(limit).map(|()| Self(molframe::fx::Comparison::AtLeast(limit)))
    }

    #[staticmethod]
    fn between(low: f64, high: f64) -> PyResult<Self> {
        finite_limit(low)?;
        finite_limit(high)?;
        if low > high {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "between requires low <= high",
            ));
        }
        Ok(Self(molframe::fx::Comparison::Between(low, high)))
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.0 {
            molframe::fx::Comparison::LessThan(_) => "less-than",
            molframe::fx::Comparison::AtMost(_) => "at-most",
            molframe::fx::Comparison::GreaterThan(_) => "greater-than",
            molframe::fx::Comparison::AtLeast(_) => "at-least",
            molframe::fx::Comparison::Between(_, _) => "between",
            _ => "unknown",
        }
    }

    #[getter]
    fn values(&self) -> Vec<f64> {
        match self.0 {
            molframe::fx::Comparison::LessThan(value)
            | molframe::fx::Comparison::AtMost(value)
            | molframe::fx::Comparison::GreaterThan(value)
            | molframe::fx::Comparison::AtLeast(value) => vec![value],
            molframe::fx::Comparison::Between(low, high) => vec![low, high],
            _ => Vec::new(),
        }
    }
}

#[pyclass(name = "MissingVerdict", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyMissingVerdict {
    Indeterminate,
    Fail,
}

impl From<PyMissingVerdict> for molframe::fx::MissingVerdict {
    fn from(value: PyMissingVerdict) -> Self {
        match value {
            PyMissingVerdict::Indeterminate => Self::Indeterminate,
            PyMissingVerdict::Fail => Self::Fail,
        }
    }
}

impl From<molframe::fx::MissingVerdict> for PyMissingVerdict {
    fn from(value: molframe::fx::MissingVerdict) -> Self {
        match value {
            molframe::fx::MissingVerdict::Indeterminate => Self::Indeterminate,
            molframe::fx::MissingVerdict::Fail => Self::Fail,
        }
    }
}

#[pyclass(name = "VerdictRule", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyVerdictRule(pub(crate) molframe::fx::VerdictRule);

#[pymethods]
impl PyVerdictRule {
    #[new]
    fn new(metric: String, comparison: &PyComparison) -> Self {
        Self(molframe::fx::VerdictRule {
            metric: metric.into_boxed_str(),
            comparison: comparison.0,
        })
    }

    #[getter]
    fn metric(&self) -> &str {
        &self.0.metric
    }

    #[getter]
    fn comparison(&self) -> PyComparison {
        PyComparison(self.0.comparison)
    }
}

#[pyclass(name = "VerdictProfile", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyVerdictProfile(pub(crate) molframe::fx::VerdictProfile);

#[pymethods]
impl PyVerdictProfile {
    #[new]
    fn new(id: String, rules: Vec<PyVerdictRule>, missing: PyMissingVerdict) -> PyResult<Self> {
        if id.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "verdict profile id must not be empty",
            ));
        }
        Ok(Self(molframe::fx::VerdictProfile::new(
            id,
            rules.into_iter().map(|rule| rule.0),
            missing.into(),
        )))
    }

    #[getter]
    fn id(&self) -> &str {
        self.0.id()
    }

    #[getter]
    fn rules(&self) -> Vec<PyVerdictRule> {
        self.0.rules().iter().cloned().map(PyVerdictRule).collect()
    }

    fn decide(&self, metrics: BTreeMap<String, f64>) -> PyVerdict {
        let metrics = metrics
            .into_iter()
            .map(|(name, value)| (name.into_boxed_str(), value))
            .collect();
        PyVerdict(self.0.decide(&metrics))
    }
}

#[pyclass(name = "RuleOutcome", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRuleOutcome {
    #[pyo3(get)]
    metric: String,
    #[pyo3(get)]
    value: Option<f64>,
    #[pyo3(get)]
    passed: Option<bool>,
}

impl From<molframe::fx::RuleOutcome> for PyRuleOutcome {
    fn from(value: molframe::fx::RuleOutcome) -> Self {
        Self {
            metric: value.metric.to_string(),
            value: value.value,
            passed: value.passed,
        }
    }
}

#[pyclass(name = "VerdictStatus", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyVerdictStatus {
    Pass,
    Fail,
    Indeterminate,
}

impl From<molframe::fx::VerdictStatus> for PyVerdictStatus {
    fn from(value: molframe::fx::VerdictStatus) -> Self {
        match value {
            molframe::fx::VerdictStatus::Pass => Self::Pass,
            molframe::fx::VerdictStatus::Fail => Self::Fail,
            molframe::fx::VerdictStatus::Indeterminate => Self::Indeterminate,
        }
    }
}

#[pyclass(name = "Verdict", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyVerdict(pub(crate) molframe::fx::Verdict);

#[pymethods]
impl PyVerdict {
    #[getter]
    fn profile(&self) -> &str {
        &self.0.profile
    }

    #[getter]
    fn status(&self) -> PyVerdictStatus {
        self.0.status.into()
    }

    #[getter]
    fn outcomes(&self) -> Vec<PyRuleOutcome> {
        self.0.outcomes.iter().cloned().map(Into::into).collect()
    }
}

fn finite_limit(value: f64) -> PyResult<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(pyo3::exceptions::PyValueError::new_err(
            "comparison thresholds must be finite",
        ))
    }
}
