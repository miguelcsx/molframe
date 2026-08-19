//! Typed bindings for the optional DSSP executable adapter.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "DsspSegment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDsspSegment {
    #[pyo3(get)]
    kind: String,
    #[pyo3(get)]
    begin_chain: String,
    #[pyo3(get)]
    begin_sequence: i64,
    #[pyo3(get)]
    end_chain: String,
    #[pyo3(get)]
    end_sequence: i64,
}

#[pyfunction]
pub(crate) fn parse_dssp_output(py: Python<'_>, output: &[u8]) -> PyResult<Vec<PyDsspSegment>> {
    py.detach(move || pdbiox::analysis::parse_dssp_output(output))
        .map(|values| values.into_iter().map(PyDsspSegment::from).collect())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn run_dssp(
    py: Python<'_>,
    executable: String,
    input: String,
) -> PyResult<Vec<PyDsspSegment>> {
    py.detach(move || pdbiox::analysis::run_dssp(executable, input))
        .map(|values| values.into_iter().map(PyDsspSegment::from).collect())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

impl From<pdbiox::analysis::DsspSegment> for PyDsspSegment {
    fn from(value: pdbiox::analysis::DsspSegment) -> Self {
        Self {
            kind: value.kind.into(),
            begin_chain: value.begin_chain.into(),
            begin_sequence: value.begin_sequence,
            end_chain: value.end_chain.into(),
            end_sequence: value.end_sequence,
        }
    }
}
