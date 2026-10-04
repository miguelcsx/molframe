//! Mechanical adapters for writing structures to mmCIF, PDB and `BinaryCIF`.

use crate::bindings::{PyStructure, findings_error};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Renders the structure as mmCIF text.
#[pyfunction]
fn to_mmcif(py: Python<'_>, structure: &PyStructure) -> PyResult<String> {
    let structure = structure.inner.clone();
    py.detach(move || molframe::write_mmcif(&structure))
        .map_err(|findings| findings_error(&findings))
}

/// Renders the structure as legacy PDB text, refusing what the format cannot hold.
#[pyfunction]
#[pyo3(signature = (structure, *, hybrid36=false, chain_map=None))]
fn to_pdb(
    py: Python<'_>,
    structure: &PyStructure,
    hybrid36: bool,
    chain_map: Option<BTreeMap<String, String>>,
) -> PyResult<String> {
    let structure = structure.inner.clone();
    let mut options = molframe::formats::pdb::PdbOptions::new().hybrid36(hybrid36);
    if let Some(chain_map) = chain_map {
        for (from, to) in chain_map {
            options = options.chain_map(from, to);
        }
    }
    py.detach(move || molframe::write_pdb(&structure, &options))
        .map_err(|findings| findings_error(&findings))
}

/// Renders the structure as deterministic `BinaryCIF` bytes.
#[pyfunction]
fn to_bcif<'py>(py: Python<'py>, structure: &PyStructure) -> PyResult<Bound<'py, PyBytes>> {
    let structure = structure.inner.clone();
    let bytes = py
        .detach(move || molframe::write_bcif(&structure))
        .map_err(|findings| findings_error(&findings))?;
    Ok(PyBytes::new(py, &bytes))
}

/// Writes the structure to `path`, choosing the format from `format` or the
/// file name.
#[pyfunction]
#[pyo3(signature = (structure, path, *, format=None))]
fn write(
    py: Python<'_>,
    structure: &PyStructure,
    path: PathBuf,
    format: Option<&str>,
) -> PyResult<()> {
    let format = match format {
        Some(name) => Some(molframe::Format::from_extension(name).ok_or_else(|| {
            crate::error::from_diagnostic(
                &molframe::Diagnostic::new(molframe::Code::E1001).with_context("name", name),
            )
        })?),
        None => None,
    };
    let structure = structure.inner.clone();
    let options = molframe::WriteOptions::canonical();
    py.detach(move || match format {
        Some(format) => molframe::write_as(&path, &structure, format, &options),
        None => molframe::write_with_options(&path, &structure, &options),
    })
    .map_err(|findings| findings_error(&findings))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(to_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(to_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(to_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)
}
