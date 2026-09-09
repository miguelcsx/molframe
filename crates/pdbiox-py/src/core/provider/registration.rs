//! Registration for native provider types.

use super::chunks::{
    PyAtomEndpoint, PyBondChunk, PyBondChunkRecord, PyFrameChunk, PyPropertyChunk, PyPropertyKind,
    PyPropertyValue, PyStructureChunk,
};
use super::metadata::{
    PyChunkDescriptor, PyChunkLayout, PyDatasetCatalog, PyDatasetDescriptor, PyPayloadKind,
};
use super::sources::{
    PyBondChunkProvider, PyFrameChunkProvider, PyPropertyChunkProvider, PyStructureChunkProvider,
};
use pyo3::prelude::*;

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("TARGET_CHUNK_BONDS", pdbiox::TARGET_CHUNK_BONDS)?;
    module.add_class::<PyPayloadKind>()?;
    module.add_class::<PyChunkLayout>()?;
    module.add_class::<PyChunkDescriptor>()?;
    module.add_class::<PyDatasetDescriptor>()?;
    module.add_class::<PyDatasetCatalog>()?;
    module.add_class::<PyPropertyKind>()?;
    module.add_class::<PyPropertyValue>()?;
    module.add_class::<PyAtomEndpoint>()?;
    module.add_class::<PyBondChunkRecord>()?;
    module.add_class::<PyBondChunk>()?;
    module.add_class::<PyStructureChunk>()?;
    module.add_class::<PyPropertyChunk>()?;
    module.add_class::<PyFrameChunk>()?;
    module.add_class::<PyStructureChunkProvider>()?;
    module.add_class::<PyBondChunkProvider>()?;
    module.add_class::<PyPropertyChunkProvider>()?;
    module.add_class::<PyFrameChunkProvider>()
}
