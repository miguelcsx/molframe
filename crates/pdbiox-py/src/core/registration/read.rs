//! Top-level structure-read entry point.

use crate::errors::read_error;
use crate::io::PyReadOptions;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyfunction]
#[pyo3(signature = (path, options=None))]
pub(crate) fn read(
    py: Python<'_>,
    path: PathBuf,
    options: Option<&PyReadOptions>,
) -> PyResult<PyStructure> {
    let options = options.map(|value| value.0.clone());
    py.detach(move || match options {
        Some(options) => pdbiox::read_with_options(path, &options).map(|value| value.0),
        None => pdbiox::read(path),
    })
    .map(PyStructure::new)
    .map_err(|findings| read_error(py, &findings))
}
