use pyo3::PyErr;
use pyo3::exceptions::PyValueError;

pub(crate) fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
