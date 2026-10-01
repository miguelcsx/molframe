//! Mechanical adapters for writing structures to mmCIF, PDB and `BinaryCIF`.

use crate::bindings::{PyStructure, findings_error};
use pyo3::types::PyBytes;
use pyo3::{exceptions::PyValueError, prelude::*};
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

#[derive(Clone, Copy)]
enum Format {
    Mmcif,
    Pdb,
    Bcif,
}

fn format_of(path: &std::path::Path, explicit: Option<&str>) -> PyResult<Format> {
    let name = match explicit {
        Some(name) => name.to_ascii_lowercase(),
        None => path
            .extension()
            .and_then(|extension| extension.to_str())
            .map(str::to_ascii_lowercase)
            .ok_or_else(|| PyValueError::new_err("cannot infer a format; pass format="))?,
    };
    match name.as_str() {
        "cif" | "mmcif" => Ok(Format::Mmcif),
        "pdb" | "ent" => Ok(Format::Pdb),
        "bcif" => Ok(Format::Bcif),
        other => Err(PyValueError::new_err(format!(
            "unsupported format {other:?}; expected mmcif, pdb or bcif"
        ))),
    }
}

/// Writes the structure to `path`, choosing the format from `format` or the
/// file extension.
#[pyfunction]
#[pyo3(signature = (structure, path, *, format=None))]
fn write(
    py: Python<'_>,
    structure: &PyStructure,
    path: PathBuf,
    format: Option<&str>,
) -> PyResult<()> {
    let format = format_of(&path, format)?;
    let structure = structure.inner.clone();
    py.detach(move || {
        let bytes = match format {
            Format::Mmcif => molframe::write_mmcif(&structure).map(String::into_bytes),
            Format::Pdb => {
                molframe::write_pdb(&structure, &molframe::formats::pdb::PdbOptions::new())
                    .map(String::into_bytes)
            }
            Format::Bcif => molframe::write_bcif(&structure),
        }
        .map_err(|findings| findings_error(&findings))?;
        std::fs::write(&path, bytes).map_err(PyErr::from)
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(to_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(to_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(to_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)
}
