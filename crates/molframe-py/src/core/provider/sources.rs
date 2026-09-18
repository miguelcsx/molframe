//! Declarative provider adapters over native structure snapshots.

use super::chunks::{PyBondChunk, PyFrameChunk, PyPropertyChunk, PyStructureChunk};
use super::metadata::{PyDatasetDescriptor, index_error, value_error};
use crate::index::{PyChunkId, PyDatasetId, PyLogicalRow, PyModelIndex};
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "BondChunkProvider", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBondChunkProvider(molframe::BondChunkProvider);

#[pymethods]
impl PyBondChunkProvider {
    #[new]
    fn new(
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        atom_dataset: PyDatasetId,
        atom_logical_start: PyLogicalRow,
        structure: &PyStructure,
    ) -> PyResult<Self> {
        molframe::BondChunkProvider::new(
            dataset.0,
            first_chunk.0,
            atom_dataset.0,
            atom_logical_start.0,
            structure.structure().clone(),
        )
        .map(Self)
        .map_err(value_error)
    }

    #[staticmethod]
    fn with_rows_per_chunk(
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        atom_dataset: PyDatasetId,
        atom_logical_start: PyLogicalRow,
        structure: &PyStructure,
        rows_per_chunk: u32,
    ) -> PyResult<Self> {
        molframe::BondChunkProvider::with_rows_per_chunk(
            dataset.0,
            first_chunk.0,
            atom_dataset.0,
            atom_logical_start.0,
            structure.structure().clone(),
            rows_per_chunk,
        )
        .map(Self)
        .map_err(value_error)
    }

    #[getter]
    fn dataset(&self) -> PyDatasetDescriptor {
        PyDatasetDescriptor(self.0.dataset())
    }

    fn chunk(&self, chunk: PyChunkId) -> PyResult<PyBondChunk> {
        self.0.chunk(chunk.0).map(PyBondChunk).map_err(index_error)
    }
}

#[pyclass(name = "StructureChunkProvider", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructureChunkProvider(molframe::StructureChunkProvider);

#[pymethods]
impl PyStructureChunkProvider {
    #[new]
    fn new(
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        structure: &PyStructure,
    ) -> PyResult<Self> {
        molframe::StructureChunkProvider::new(
            dataset.0,
            first_chunk.0,
            structure.structure().clone(),
        )
        .map(Self)
        .map_err(value_error)
    }

    #[getter]
    fn dataset(&self) -> PyDatasetDescriptor {
        PyDatasetDescriptor(self.0.dataset())
    }

    fn chunk(&self, chunk: PyChunkId) -> PyResult<PyStructureChunk> {
        self.0
            .chunk(chunk.0)
            .map(PyStructureChunk)
            .map_err(index_error)
    }
}

#[pyclass(name = "PropertyChunkProvider", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPropertyChunkProvider(molframe::PropertyChunkProvider);

#[pymethods]
impl PyPropertyChunkProvider {
    #[staticmethod]
    fn plddt(
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        structure: &PyStructure,
    ) -> PyResult<Self> {
        molframe::PropertyChunkProvider::plddt(
            dataset.0,
            first_chunk.0,
            structure.structure().clone(),
        )
        .map(Self)
        .map_err(value_error)
    }

    #[new]
    fn new(
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        structure: &PyStructure,
        name: String,
    ) -> PyResult<Self> {
        molframe::PropertyChunkProvider::new(
            dataset.0,
            first_chunk.0,
            structure.structure().clone(),
            name,
        )
        .map(Self)
        .map_err(value_error)
    }

    #[getter]
    fn dataset(&self) -> PyDatasetDescriptor {
        PyDatasetDescriptor(self.0.dataset())
    }

    fn chunk(&self, chunk: PyChunkId) -> PyResult<PyPropertyChunk> {
        self.0
            .chunk(chunk.0)
            .map(PyPropertyChunk)
            .map_err(index_error)
    }
}

#[pyclass(name = "FrameChunkProvider", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFrameChunkProvider(molframe::FrameChunkProvider);

#[pymethods]
impl PyFrameChunkProvider {
    #[new]
    fn new(
        dataset: PyDatasetId,
        first_chunk: PyChunkId,
        structure: &PyStructure,
        model: PyModelIndex,
    ) -> PyResult<Self> {
        molframe::FrameChunkProvider::new(
            dataset.0,
            first_chunk.0,
            structure.structure().clone(),
            model.0,
        )
        .map(Self)
        .map_err(value_error)
    }

    #[getter]
    fn dataset(&self) -> PyDatasetDescriptor {
        PyDatasetDescriptor(self.0.dataset())
    }

    fn chunk(&self, chunk: PyChunkId) -> PyResult<PyFrameChunk> {
        self.0.chunk(chunk.0).map(PyFrameChunk).map_err(index_error)
    }
}
