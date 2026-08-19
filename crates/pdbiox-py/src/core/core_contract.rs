//! Typed projections of the core reproducibility and comparison contracts.

use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::BTreeMap;
use std::path::PathBuf;

type NativeParameter = pdbiox::core::contract::ParameterValue;

#[pyclass(name = "AlgorithmId", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAlgorithmId(pub(crate) pdbiox::AlgorithmId);

#[pymethods]
impl PyAlgorithmId {
    #[new]
    fn new(name: String, version: String) -> Self {
        Self(pdbiox::AlgorithmId::new(name, version))
    }

    #[getter]
    fn name(&self) -> &str {
        self.0.name()
    }

    #[getter]
    fn version(&self) -> &str {
        self.0.version()
    }

    fn __repr__(&self) -> String {
        format!("AlgorithmId({:?}, {:?})", self.name(), self.version())
    }
}

#[pyclass(name = "DictionaryVersion", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDictionaryVersion(pub(crate) pdbiox::DictionaryVersion);

#[pymethods]
impl PyDictionaryVersion {
    #[new]
    fn new(version: String) -> Self {
        Self(pdbiox::DictionaryVersion::new(version))
    }

    #[getter]
    fn value(&self) -> &str {
        self.0.as_str()
    }

    fn __str__(&self) -> &str {
        self.0.as_str()
    }
}

#[pyclass(name = "ProfileId", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyProfileId(pub(crate) pdbiox::ProfileId);

#[pymethods]
impl PyProfileId {
    #[staticmethod]
    fn default() -> Self {
        Self(pdbiox::ProfileId::DEFAULT)
    }

    #[getter]
    fn value(&self) -> &str {
        self.0.as_str()
    }

    fn __str__(&self) -> &str {
        self.0.as_str()
    }
}

#[pyclass(name = "Fingerprint", frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyFingerprint(pub(crate) pdbiox::core::contract::Fingerprint);

#[pymethods]
impl PyFingerprint {
    #[new]
    fn new(data: &[u8]) -> Self {
        Self(pdbiox::core::contract::Fingerprint::of(data))
    }

    #[staticmethod]
    fn of(data: &[u8]) -> Self {
        Self(pdbiox::core::contract::Fingerprint::of(data))
    }

    #[getter]
    fn value(&self) -> u64 {
        self.0.get()
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }
}

#[pyclass(name = "ParameterValue", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyParameterValue(pub(crate) NativeParameter);

#[pymethods]
impl PyParameterValue {
    #[staticmethod]
    fn boolean(value: bool) -> Self {
        Self(NativeParameter::Boolean(value))
    }

    #[staticmethod]
    fn integer(value: i64) -> Self {
        Self(NativeParameter::Integer(value))
    }

    #[staticmethod]
    fn float(value: f64) -> PyResult<Self> {
        NativeParameter::finite_float(value)
            .map(Self)
            .ok_or_else(|| PyValueError::new_err("parameter float must be finite"))
    }

    #[staticmethod]
    fn text(value: String) -> Self {
        Self(NativeParameter::Text(value.into()))
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.0 {
            NativeParameter::Boolean(_) => "boolean",
            NativeParameter::Integer(_) => "integer",
            NativeParameter::Float(_) => "float",
            NativeParameter::Text(_) => "text",
        }
    }

    #[getter]
    fn bool_value(&self) -> Option<bool> {
        match self.0 {
            NativeParameter::Boolean(value) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn integer_value(&self) -> Option<i64> {
        match self.0 {
            NativeParameter::Integer(value) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn float_value(&self) -> Option<f64> {
        match self.0 {
            NativeParameter::Float(bits) => Some(f64::from_bits(bits)),
            _ => None,
        }
    }

    #[getter]
    fn text_value(&self) -> Option<String> {
        match &self.0 {
            NativeParameter::Text(value) => Some(value.to_string()),
            _ => None,
        }
    }

    fn to_python(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match &self.0 {
            NativeParameter::Boolean(value) => {
                Ok(value.into_pyobject(py)?.to_owned().into_any().unbind())
            }
            NativeParameter::Integer(value) => Ok(value.into_pyobject(py)?.into_any().unbind()),
            NativeParameter::Float(bits) => {
                let value = infallible(f64::from_bits(*bits).into_pyobject(py));
                Ok(value.into_any().unbind())
            }
            NativeParameter::Text(value) => {
                let value = infallible(value.to_string().into_pyobject(py));
                Ok(value.into_any().unbind())
            }
        }
    }
}

fn infallible<T>(value: Result<T, std::convert::Infallible>) -> T {
    match value {
        Ok(value) => value,
        Err(error) => match error {},
    }
}

#[pyclass(name = "AnalysisParameters", from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyAnalysisParameters(pub(crate) pdbiox::AnalysisParameters);

#[pymethods]
impl PyAnalysisParameters {
    #[new]
    #[pyo3(signature = (values=None))]
    fn new(values: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut parameters = BTreeMap::new();
        if let Some(values) = values {
            for (key, value) in values.iter() {
                let key = key.extract::<String>()?;
                let value = value.extract::<PyParameterValue>()?;
                parameters.insert(key.into_boxed_str(), value.0);
            }
        }
        Ok(Self(parameters))
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __contains__(&self, name: &str) -> bool {
        self.0.contains_key(name)
    }

    fn get(&self, name: &str) -> Option<PyParameterValue> {
        self.0.get(name).cloned().map(PyParameterValue)
    }

    fn __getitem__(&self, name: &str) -> PyResult<PyParameterValue> {
        self.get(name)
            .ok_or_else(|| PyKeyError::new_err(name.to_owned()))
    }

    fn set(&mut self, name: String, value: PyParameterValue) {
        self.0.insert(name.into_boxed_str(), value.0);
    }

    fn items(&self) -> Vec<(String, PyParameterValue)> {
        self.0
            .iter()
            .map(|(name, value)| (name.to_string(), PyParameterValue(value.clone())))
            .collect()
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        for (name, value) in &self.0 {
            result.set_item(
                name.to_string(),
                PyParameterValue(value.clone()).to_python(py)?,
            )?;
        }
        Ok(result)
    }
}

#[pyclass(name = "SourceRef", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySourceRef(pub(crate) pdbiox::SourceRef);

#[pymethods]
impl PySourceRef {
    #[staticmethod]
    fn none() -> Self {
        Self(pdbiox::SourceRef::None)
    }

    #[staticmethod]
    fn path(path: PathBuf) -> Self {
        Self(pdbiox::SourceRef::path(path))
    }

    #[staticmethod]
    fn url(url: String) -> Self {
        Self(pdbiox::SourceRef::Url(url.into_boxed_str()))
    }

    #[staticmethod]
    fn memory() -> Self {
        Self(pdbiox::SourceRef::Memory)
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match self.0 {
            pdbiox::SourceRef::None => "none",
            pdbiox::SourceRef::Path(_) => "path",
            pdbiox::SourceRef::Url(_) => "url",
            pdbiox::SourceRef::Memory => "memory",
            _ => "unknown",
        }
    }

    #[getter]
    fn path_value(&self) -> Option<String> {
        self.0
            .as_path()
            .map(|path| path.to_string_lossy().into_owned())
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }
}

#[pyclass(name = "StructureDifferenceOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyStructureDifferenceOptions(pub(crate) pdbiox::StructureDifferenceOptions);

#[pymethods]
impl PyStructureDifferenceOptions {
    #[new]
    #[pyo3(signature = (coordinate_tolerance=0.0))]
    fn new(coordinate_tolerance: f32) -> Self {
        Self(pdbiox::StructureDifferenceOptions {
            coordinate_tolerance,
        })
    }

    #[getter]
    fn coordinate_tolerance(&self) -> f32 {
        self.0.coordinate_tolerance
    }
}

#[pyclass(name = "ValueDifference", frozen, skip_from_py_object)]
pub(crate) struct PyValueDifference {
    left: Py<PyAny>,
    right: Py<PyAny>,
}

#[pymethods]
impl PyValueDifference {
    #[new]
    fn new(left: Py<PyAny>, right: Py<PyAny>) -> Self {
        Self { left, right }
    }

    #[getter]
    fn left(&self, py: Python<'_>) -> Py<PyAny> {
        self.left.clone_ref(py)
    }

    #[getter]
    fn right(&self, py: Python<'_>) -> Py<PyAny> {
        self.right.clone_ref(py)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAlgorithmId>()?;
    module.add_class::<PyDictionaryVersion>()?;
    module.add_class::<PyProfileId>()?;
    module.add_class::<PyFingerprint>()?;
    module.add_class::<PyParameterValue>()?;
    module.add_class::<PyAnalysisParameters>()?;
    module.add_class::<PySourceRef>()?;
    module.add_class::<PyStructureDifferenceOptions>()?;
    module.add_class::<PyValueDifference>()?;
    Ok(())
}
