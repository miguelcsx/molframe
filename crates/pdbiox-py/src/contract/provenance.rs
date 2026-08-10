//! Python projection of diagnostics and reproducibility metadata.

use pdbiox::{Diagnostic, ParameterValue, Provenance};
use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "Diagnostic", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDiagnostic {
    #[pyo3(get)]
    code: String,
    #[pyo3(get)]
    message: String,
    #[pyo3(get)]
    remedy: String,
    #[pyo3(get)]
    span: Option<(u32, u32, u32, u32)>,
}

impl From<Diagnostic> for PyDiagnostic {
    fn from(value: Diagnostic) -> Self {
        Self {
            code: value.code().to_string(),
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
