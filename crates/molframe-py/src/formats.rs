//! Mechanical adapters for writing structures to mmCIF, PDB and `BinaryCIF`.

use crate::bindings::{PyStructure, findings_error};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// The mmCIF writing decisions a call states.
///
/// Connection identifiers are generated unless the caller turns them off, so a
/// structure with bonds writes by default; the block name falls back to the
/// entry id and then to a fixed name inside the writer.
fn cif_options(
    block_id: Option<&str>,
    generated_connections: bool,
) -> molframe::formats::cif::CifWriteOptions {
    let mut options = molframe::formats::cif::CifWriteOptions::new();
    if let Some(block_id) = block_id {
        options = options.with_block_id(block_id);
    }
    if !generated_connections {
        options = options.without_generated_connection_ids();
    }
    options
}

/// Renders the structure as mmCIF text.
///
/// The space group, resolution, assemblies, non-crystallographic operators,
/// secondary structure and connections the structure holds are written.
#[pyfunction]
#[pyo3(signature = (structure, *, block_id=None, generated_connections=true))]
fn to_mmcif(
    py: Python<'_>,
    structure: &PyStructure,
    block_id: Option<&str>,
    generated_connections: bool,
) -> PyResult<String> {
    let structure = structure.inner.clone();
    let options = cif_options(block_id, generated_connections);
    py.detach(move || molframe::write_mmcif_with_options(&structure, &options))
        .map_err(|findings| findings_error(&findings))
}

/// Copies an mmCIF document to `destination` category for category, in order.
///
/// This is the lossless route: nothing is lowered to a structure, so nothing the
/// structure model does not hold can be dropped.
#[pyfunction]
fn preserve_mmcif(py: Python<'_>, source: PathBuf, destination: PathBuf) -> PyResult<()> {
    let text = py
        .detach(move || {
            molframe::read_document(&source)
                .map(|document| molframe::formats::cif::write_preserving(&document))
        })
        .map_err(|findings| findings_error(&findings))?;
    std::fs::write(&destination, text).map_err(move |error| {
        crate::error::kernel(
            molframe::Diagnostic::new(molframe::Code::E7101)
                .with_message(format!("{}: {error}", destination.display())),
        )
    })
}

/// Renders the structure as legacy PDB text, refusing what the format cannot hold.
#[pyfunction]
#[pyo3(signature = (structure, *, hybrid36=false, chain_map=None, namespace="auto"))]
fn to_pdb(
    py: Python<'_>,
    structure: &PyStructure,
    hybrid36: bool,
    chain_map: Option<BTreeMap<String, String>>,
    namespace: &str,
) -> PyResult<String> {
    let options = pdb_options(&structure.inner, hybrid36, chain_map, namespace)?;
    let structure = structure.inner.clone();
    py.detach(move || molframe::write_pdb(&structure, &options))
        .map_err(|findings| findings_error(&findings))
}

/// Renders the structure as deterministic `BinaryCIF` bytes.
#[pyfunction]
#[pyo3(signature = (structure, *, block_id=None, generated_connections=true, transport=false))]
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
    let options = cif_options(block_id, generated_connections);
    let bytes = py
        .detach(move || molframe::write_bcif_with_options(&structure, &options))
        .map_err(|findings| findings_error(&findings))?;
    Ok(PyBytes::new(py, &bytes))
}

/// The PDB writing decisions a call states: hybrid-36 numbering, chain renames
/// and the identifier namespace.
///
/// `auto` writes the author identifiers the format was designed around when
/// every chain and residue has them, and the normalised label identifiers
/// otherwise. Label sequence identifiers are absent for ligands and water, so
/// choosing them by default made most deposited files unwritable.
fn pdb_options(
    structure: &molframe::Structure,
    hybrid36: bool,
    chain_map: Option<BTreeMap<String, String>>,
    namespace: &str,
) -> PyResult<molframe::formats::pdb::PdbOptions> {
    use molframe::formats::pdb::PdbIdentifierNamespace as Namespace;
    let namespace = match namespace {
        "auto" => {
            let data = structure.engine().data();
            let has_author = data.chains().all(|chain| chain.auth_label().is_some())
                && data
                    .residues()
                    .all(|residue| residue.auth_seq_id().is_some());
            if has_author {
                Namespace::Auth
            } else {
                Namespace::Label
            }
        }
        "auth" => Namespace::Auth,
        "label" => Namespace::Label,
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "namespace must be 'auto', 'auth' or 'label', not {other:?}"
            )));
        }
    };
    let mut options = molframe::formats::pdb::PdbOptions::new()
        .hybrid36(hybrid36)
        .namespace(namespace);
    for (from, to) in chain_map.into_iter().flatten() {
        options = options.chain_map(from, to);
    }
    Ok(options)
}

/// Writes the structure to `path`, choosing the format from `format` or the
/// file name.
///
/// A `.gz` or `.zst` suffix applies deterministic compression. `hybrid36` and
/// `chain_map` and `namespace` are the PDB decisions of `to_pdb`; `block_id` and
/// `generated_connections` are the `mmCIF` and `BinaryCIF` decisions of `to_mmcif`;
/// `memory_limit` bounds, in bytes, the working memory the writer may keep.
#[pyfunction]
#[pyo3(signature = (structure, path, *, format=None, hybrid36=false, chain_map=None, namespace="auto", block_id=None, generated_connections=true, memory_limit=None))]
fn write(
    py: Python<'_>,
    structure: &PyStructure,
    path: PathBuf,
    format: Option<&str>,
    hybrid36: bool,
    chain_map: Option<BTreeMap<String, String>>,
    namespace: &str,
    block_id: Option<&str>,
    generated_connections: bool,
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
    let mut options = molframe::WriteOptions::canonical()
        .with_pdb(pdb_options(
            &structure.inner,
            hybrid36,
            chain_map,
            namespace,
        )?)
        .with_cif(cif_options(block_id, generated_connections));
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
    module.add_function(wrap_pyfunction!(preserve_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(to_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(to_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)
}
