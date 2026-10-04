//! Interoperability: Arrow streams, `DLPack` coordinates, graph tensors and datasets.

mod dataset;
mod dlpack;
mod graph;
mod stream;
mod tables;

pub(crate) use stream::arrow_capsule;
pub(crate) use tables::{PyBondTable, columns_stream};

use crate::bindings::PyStructure;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Writes the atom table as a Parquet file, one record batch at a time.
///
/// `metadata` is stored in the file's schema.
#[pyfunction]
#[pyo3(signature = (structure, path, *, metadata=BTreeMap::new()))]
fn write_atoms_parquet(
    py: Python<'_>,
    structure: &PyStructure,
    path: PathBuf,
    metadata: BTreeMap<String, String>,
) -> PyResult<()> {
    let source = structure.inner.clone();
    py.detach(move || {
        molframe::interop::write_atom_parquet_with_metadata(&path, source.engine(), metadata)
    })
    .map_err(crate::error::kernel)
}

/// Writes the atom table as an Arrow IPC file, one record batch at a time.
///
/// `metadata` is stored in the file's schema.
#[pyfunction]
#[pyo3(signature = (structure, path, *, metadata=BTreeMap::new()))]
fn write_atoms_ipc(
    py: Python<'_>,
    structure: &PyStructure,
    path: PathBuf,
    metadata: BTreeMap<String, String>,
) -> PyResult<()> {
    let source = structure.inner.clone();
    py.detach(move || {
        molframe::interop::write_atom_ipc_with_metadata(&path, source.engine(), metadata)
    })
    .map_err(crate::error::kernel)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<dataset::PyDataset>()?;
    module.add_class::<dataset::PyDatasetSplit>()?;
    module.add_class::<dataset::PyManifestEntry>()?;
    module.add_class::<dlpack::PyCoordinateTensor>()?;
    module.add_class::<graph::PyGraph>()?;
    module.add_function(wrap_pyfunction!(dlpack::coordinates, module)?)?;
    module.add_function(wrap_pyfunction!(graph::graph, module)?)?;
    module.add_function(wrap_pyfunction!(write_atoms_ipc, module)?)?;
    module.add_function(wrap_pyfunction!(write_atoms_parquet, module)?)
}
