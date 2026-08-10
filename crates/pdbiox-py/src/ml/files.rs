//! Native Arrow and Parquet atom-table writers.

use crate::structure::PyStructure;
use pyo3::exceptions::PyOSError;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyfunction]
pub(crate) fn write_atom_ipc(
    py: Python<'_>,
    path: PathBuf,
    structure: &PyStructure,
) -> PyResult<()> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_atom_ipc(path, &structure))
        .map_err(io_error)
}

#[pyfunction]
pub(crate) fn write_atom_parquet(
    py: Python<'_>,
    path: PathBuf,
    structure: &PyStructure,
) -> PyResult<()> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_atom_parquet(path, &structure))
        .map_err(io_error)
}

fn io_error(error: impl std::fmt::Display) -> PyErr {
    PyOSError::new_err(error.to_string())
}
