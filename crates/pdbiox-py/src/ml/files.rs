//! Native Arrow and Parquet atom-table writers.

use crate::structure::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[pyfunction]
pub(crate) fn write_atom_ipc(
    py: Python<'_>,
    path: PathBuf,
    structure: &PyStructure,
) -> PyResult<()> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_atom_ipc(path, &structure))
        .map_err(|error| crate::errors::table_file_error(&error))
}

#[pyfunction]
pub(crate) fn write_atom_parquet(
    py: Python<'_>,
    path: PathBuf,
    structure: &PyStructure,
) -> PyResult<()> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_atom_parquet(path, &structure))
        .map_err(|error| crate::errors::table_file_error(&error))
}

#[pyfunction]
#[pyo3(signature = (path, structure, *, metadata=None))]
pub(crate) fn write_atom_ipc_with_metadata(
    py: Python<'_>,
    path: PathBuf,
    structure: &PyStructure,
    metadata: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let metadata = metadata_map(metadata)?;
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_atom_ipc_with_metadata(path, &structure, metadata))
        .map_err(|error| crate::errors::table_file_error(&error))
}

#[pyfunction]
#[pyo3(signature = (path, structure, *, metadata=None))]
pub(crate) fn write_atom_parquet_with_metadata(
    py: Python<'_>,
    path: PathBuf,
    structure: &PyStructure,
    metadata: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let metadata = metadata_map(metadata)?;
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_atom_parquet_with_metadata(path, &structure, metadata))
        .map_err(|error| crate::errors::table_file_error(&error))
}

fn metadata_map(metadata: Option<&Bound<'_, PyDict>>) -> PyResult<BTreeMap<String, String>> {
    let Some(metadata) = metadata else {
        return Ok(BTreeMap::new());
    };
    metadata
        .iter()
        .map(|(key, value)| Ok((key.extract::<String>()?, value.extract::<String>()?)))
        .collect()
}
