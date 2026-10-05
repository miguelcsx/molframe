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
    let options = pdb_options(hybrid36, chain_map);
    py.detach(move || molframe::write_pdb(&structure, &options))
        .map_err(|findings| findings_error(&findings))
}

/// Renders the structure as deterministic `BinaryCIF` bytes.
#[pyfunction]
#[pyo3(signature = (structure, *, block_id=None, generated_connections=false, transport=false))]
fn to_bcif<'py>(
    py: Python<'py>,
    structure: &PyStructure,
    block_id: Option<&str>,
    generated_connections: bool,
    transport: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    let structure = if transport {
        structure
            .inner
            .with_transport_identifiers()
            .map_err(crate::error::kernel)?
    } else {
        structure.inner.clone()
    };
    let mut options = molframe::formats::cif::CifWriteOptions::new();
    if let Some(block_id) = block_id {
        options = options.with_block_id(block_id);
    }
    if generated_connections {
        options = options
            .with_generated_connection_ids()
            .with_connection_type_id("covale");
    }
    let bytes = py
        .detach(move || molframe::write_bcif_with_options(&structure, &options))
        .map_err(|findings| findings_error(&findings))?;
    Ok(PyBytes::new(py, &bytes))
}

/// The PDB writing decisions a call states: hybrid-36 numbering and chain renames.
fn pdb_options(
    hybrid36: bool,
    chain_map: Option<BTreeMap<String, String>>,
) -> molframe::formats::pdb::PdbOptions {
    let mut options = molframe::formats::pdb::PdbOptions::new().hybrid36(hybrid36);
    for (from, to) in chain_map.into_iter().flatten() {
        options = options.chain_map(from, to);
    }
    options
}

/// Writes the structure to `path`, choosing the format from `format` or the
/// file name.
///
/// A `.gz` or `.zst` suffix applies deterministic compression. `hybrid36` and
/// `chain_map` are the PDB decisions of `to_pdb`; `memory_limit` bounds, in
/// bytes, the working memory the writer may keep.
#[pyfunction]
#[pyo3(signature = (structure, path, *, format=None, hybrid36=false, chain_map=None, memory_limit=None))]
fn write(
    py: Python<'_>,
    structure: &PyStructure,
    path: PathBuf,
    format: Option<&str>,
    hybrid36: bool,
    chain_map: Option<BTreeMap<String, String>>,
    memory_limit: Option<usize>,
) -> PyResult<()> {
    let format = match format {
        Some(name) => match name
            .parse::<molframe::Format>()
            .map_err(crate::error::kernel)?
        {
            molframe::Format::Auto => None,
            chosen => Some(chosen),
        },
        None => None,
    };
    let mut options =
        molframe::WriteOptions::canonical().with_pdb(pdb_options(hybrid36, chain_map));
    if let Some(bytes) = memory_limit {
        options = options.with_output(molframe::OutputOptions::default().with_memory_limit(bytes));
    }
    let structure = structure.inner.clone();
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
