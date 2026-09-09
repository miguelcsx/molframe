//! Python projection of diagnostics and reproducibility metadata.

use super::super::core_diagnostic::{PyCode, PyContextItem, PyRendered, PySeverity, PyStrictness};
use crate::core_values::PyByteSpan;
use pdbiox::{Diagnostic, ParameterValue, Provenance};
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "Diagnostic", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDiagnostic {
    pub(crate) inner: Diagnostic,
    #[pyo3(get)]
    code: String,
    #[pyo3(get)]
    message: String,
    #[pyo3(get)]
    remedy: String,
    #[pyo3(get)]
    span: Option<(u64, u64, u64, u64)>,
}

impl From<Diagnostic> for PyDiagnostic {
    fn from(value: Diagnostic) -> Self {
        let code = value.code();
        Self {
            inner: value.clone(),
            code: code.to_string(),
            message: value.message().to_owned(),
            remedy: value.remedy().to_owned(),
            span: value.span().map(|span| {
                (
                    span.start.byte_offset,
                    span.end,
                    span.start.line,
                    span.start.column,
                )
            }),
        }
    }
}

#[pymethods]
impl PyDiagnostic {
    #[new]
    fn new(code: &PyCode) -> Self {
        pdbiox::Diagnostic::new(code.0).into()
    }

    fn with_message(&self, message: &str) -> Self {
        let mut value = self.clone();
        value.inner = value.inner.with_message(message);
        value.refresh()
    }

    fn with_severity(&self, severity: PySeverity) -> Self {
        let mut value = self.clone();
        value.inner = value.inner.with_severity(severity.into());
        value.refresh()
    }

    fn at(&self, span: PyByteSpan) -> Self {
        let mut value = self.clone();
        value.inner = value.inner.at(span.0);
        value.refresh()
    }

    fn in_category(&self, category: &str) -> Self {
        let mut value = self.clone();
        value.inner = value.inner.in_category(category);
        value.refresh()
    }

    fn about_field(&self, field: &str) -> Self {
        let mut value = self.clone();
        value.inner = value.inner.about_field(field);
        value.refresh()
    }

    fn at_row(&self, row: u64) -> Self {
        let mut value = self.clone();
        value.inner = value.inner.at_row(row);
        value.refresh()
    }

    #[getter]
    fn code_value(&self) -> PyCode {
        PyCode(self.inner.code())
    }

    #[getter]
    fn severity_value(&self) -> PySeverity {
        self.inner.severity().into()
    }

    #[getter]
    fn category_value(&self) -> Option<String> {
        self.inner.category().map(str::to_owned)
    }

    #[getter]
    fn field_value(&self) -> Option<String> {
        self.inner.field().map(str::to_owned)
    }

    #[getter]
    fn row(&self) -> Option<u64> {
        self.inner.row()
    }

    #[getter]
    fn context(&self) -> Vec<PyContextItem> {
        self.inner
            .context()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn is_error(&self, strictness: PyStrictness) -> bool {
        self.inner.is_error(strictness.into())
    }

    #[pyo3(signature = (source=None, *, origin=None, color=false))]
    fn render(&self, source: Option<&[u8]>, origin: Option<&str>, color: bool) -> PyRendered {
        PyRendered::from_diagnostic(&self.inner, source, origin, color)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Diagnostic('{}')", self.inner.code())
    }
}

impl PyDiagnostic {
    fn refresh(mut self) -> Self {
        self.code = self.inner.code().to_string();
        self.message = self.inner.message().to_owned();
        self.remedy = self.inner.remedy().to_owned();
        self.span = self.inner.span().map(|span| {
            (
                span.start.byte_offset,
                span.end,
                span.start.line,
                span.start.column,
            )
        });
        self
    }
}

#[derive(Clone, Debug)]
enum PyParameter {
    Boolean(bool),
    Integer(i64),
    Float(f64),
    Text(String),
}

impl From<ParameterValue> for PyParameter {
    fn from(value: ParameterValue) -> Self {
        match value {
            ParameterValue::Boolean(value) => Self::Boolean(value),
            ParameterValue::Integer(value) => Self::Integer(value),
            ParameterValue::Float(bits) => Self::Float(f64::from_bits(bits)),
            ParameterValue::Text(value) => Self::Text(value.into()),
        }
    }
}

#[pyclass(name = "Provenance", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyProvenance {
    pub(crate) inner: Provenance,
    #[pyo3(get)]
    pdbiox_version: String,
    #[pyo3(get)]
    input_source: String,
    #[pyo3(get)]
    input_fingerprint: Option<String>,
    #[pyo3(get)]
    policy_fingerprint: String,
    #[pyo3(get)]
    profile: Option<String>,
    #[pyo3(get)]
    algorithm_name: Option<String>,
    #[pyo3(get)]
    algorithm_version: Option<String>,
    parameters: Vec<(String, PyParameter)>,
    #[pyo3(get)]
    schema_version: Option<String>,
    #[pyo3(get)]
    component_version: Option<String>,
    #[pyo3(get)]
    timestamp: Option<String>,
    #[pyo3(get)]
    fingerprint: String,
}

#[pymethods]
impl PyProvenance {
    #[getter]
    fn parameters(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let values = PyDict::new(py);
        for (name, value) in &self.parameters {
            match value {
                PyParameter::Boolean(value) => values.set_item(name, value)?,
                PyParameter::Integer(value) => values.set_item(name, value)?,
                PyParameter::Float(value) => values.set_item(name, value)?,
                PyParameter::Text(value) => values.set_item(name, value)?,
            }
        }
        Ok(values.unbind())
    }
}

impl From<Provenance> for PyProvenance {
    fn from(value: Provenance) -> Self {
        let fingerprint = value.fingerprint().to_string();
        let (algorithm_name, algorithm_version) =
            value.algorithm.as_ref().map_or((None, None), |algorithm| {
                (
                    Some(algorithm.name().to_owned()),
                    Some(algorithm.version().to_owned()),
                )
            });
        Self {
            inner: value.clone(),
            pdbiox_version: value.pdbiox_version.to_owned(),
            input_source: value.input_source.to_string(),
            input_fingerprint: value.input_fingerprint.map(|item| item.to_string()),
            policy_fingerprint: value.policy_fingerprint.to_string(),
            profile: value.profile.map(|item| item.as_str().to_owned()),
            algorithm_name,
            algorithm_version,
            parameters: value
                .parameters
                .into_iter()
                .map(|(name, parameter)| (name.into(), parameter.into()))
                .collect(),
            schema_version: value.schema_version.map(|item| item.as_str().to_owned()),
            component_version: value.component_version.map(|item| item.as_str().to_owned()),
            timestamp: value.timestamp.map(Into::into),
            fingerprint,
        }
    }
}
