//! Dataset and chunk metadata wrappers.

use crate::index::{PyChunkId, PyDatasetId, PyLocalRow, PyLogicalRow};
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

#[cfg(test)]
#[path = "metadata_tests.rs"]
mod tests;

#[pyclass(name = "PayloadKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPayloadKind {
    Structure,
    BondTopology,
    Property,
    Frame,
}

impl From<pdbiox::PayloadKind> for PyPayloadKind {
    fn from(value: pdbiox::PayloadKind) -> Self {
        match value {
            pdbiox::PayloadKind::Structure => Self::Structure,
            pdbiox::PayloadKind::BondTopology => Self::BondTopology,
            pdbiox::PayloadKind::Property => Self::Property,
            pdbiox::PayloadKind::Frame => Self::Frame,
        }
    }
}

impl From<PyPayloadKind> for pdbiox::PayloadKind {
    fn from(value: PyPayloadKind) -> Self {
        match value {
            PyPayloadKind::Structure => Self::Structure,
            PyPayloadKind::BondTopology => Self::BondTopology,
            PyPayloadKind::Property => Self::Property,
            PyPayloadKind::Frame => Self::Frame,
        }
    }
}

#[pyclass(name = "ChunkLayout", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyChunkLayout(pub(crate) pdbiox::ChunkLayout);

#[pymethods]
impl PyChunkLayout {
    #[staticmethod]
    fn regular(rows_per_chunk: u32) -> Self {
        Self(pdbiox::ChunkLayout::Regular { rows_per_chunk })
    }

    #[staticmethod]
    fn source_defined(target_rows: u32) -> Self {
        Self(pdbiox::ChunkLayout::SourceDefined { target_rows })
    }

    #[getter]
    fn target_rows(&self) -> u32 {
        self.0.target_rows()
    }

    #[getter]
    fn is_regular(&self) -> bool {
        matches!(self.0, pdbiox::ChunkLayout::Regular { .. })
    }
}

#[pyclass(name = "ChunkDescriptor", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyChunkDescriptor(pub(crate) pdbiox::ChunkDescriptor);

#[pymethods]
impl PyChunkDescriptor {
    #[new]
    fn new(
        dataset: PyDatasetId,
        chunk: PyChunkId,
        logical_start: PyLogicalRow,
        rows: u32,
    ) -> PyResult<Self> {
        pdbiox::ChunkDescriptor::new(dataset.0, chunk.0, logical_start.0, rows)
            .map(Self)
            .map_err(value_error)
    }

    #[getter]
    fn dataset(&self) -> PyDatasetId {
        PyDatasetId(self.0.dataset())
    }

    #[getter]
    fn chunk(&self) -> PyChunkId {
        PyChunkId(self.0.chunk())
    }

    #[getter]
    fn logical_start(&self) -> PyLogicalRow {
        PyLogicalRow(self.0.logical_start())
    }

    #[getter]
    fn rows(&self) -> u32 {
        self.0.rows()
    }

    fn resolve(&self, local: PyLocalRow) -> PyResult<PyLogicalRow> {
        self.0
            .resolve(local.0)
            .map(PyLogicalRow)
            .map_err(index_error)
    }
}

#[pyclass(name = "DatasetDescriptor", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDatasetDescriptor(pub(crate) pdbiox::DatasetDescriptor);

#[pymethods]
impl PyDatasetDescriptor {
    #[staticmethod]
    fn regular(
        dataset: PyDatasetId,
        payload: PyPayloadKind,
        logical_rows: u64,
        rows_per_chunk: u32,
        first_chunk: PyChunkId,
    ) -> PyResult<Self> {
        pdbiox::DatasetDescriptor::regular(
            dataset.0,
            payload.into(),
            logical_rows,
            rows_per_chunk,
            first_chunk.0,
        )
        .map(Self)
        .map_err(value_error)
    }

    #[staticmethod]
    fn source_defined(
        dataset: PyDatasetId,
        payload: PyPayloadKind,
        logical_rows: u64,
        chunk_count: u64,
        first_chunk: PyChunkId,
        target_rows: u32,
    ) -> PyResult<Self> {
        pdbiox::DatasetDescriptor::source_defined(
            dataset.0,
            payload.into(),
            logical_rows,
            chunk_count,
            first_chunk.0,
            target_rows,
        )
        .map(Self)
        .map_err(value_error)
    }

    #[getter]
    fn id(&self) -> PyDatasetId {
        PyDatasetId(self.0.id())
    }

    #[getter]
    fn payload(&self) -> PyPayloadKind {
        self.0.payload().into()
    }

    #[getter]
    fn logical_rows(&self) -> u64 {
        self.0.logical_rows()
    }

    #[getter]
    fn chunk_count(&self) -> u64 {
        self.0.chunk_count()
    }

    #[getter]
    fn first_chunk(&self) -> PyChunkId {
        PyChunkId(self.0.first_chunk())
    }

    #[getter]
    fn layout(&self) -> PyChunkLayout {
        PyChunkLayout(self.0.layout())
    }

    fn regular_chunk(&self, chunk: PyChunkId) -> PyResult<PyChunkDescriptor> {
        self.0
            .regular_chunk(chunk.0)
            .map(PyChunkDescriptor)
            .map_err(index_error)
    }
}

#[pyclass(name = "DatasetCatalog", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDatasetCatalog(pub(crate) pdbiox::DatasetCatalog);

#[pymethods]
impl PyDatasetCatalog {
    #[new]
    fn new(datasets: Vec<PyDatasetDescriptor>) -> PyResult<Self> {
        pdbiox::DatasetCatalog::new(datasets.into_iter().map(|item| item.0).collect())
            .map(Self)
            .map_err(value_error)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[getter]
    fn datasets(&self) -> Vec<PyDatasetDescriptor> {
        self.0
            .datasets()
            .iter()
            .copied()
            .map(PyDatasetDescriptor)
            .collect()
    }

    fn get(&self, dataset: PyDatasetId) -> Option<PyDatasetDescriptor> {
        self.0.get(dataset.0).map(PyDatasetDescriptor)
    }
}

pub(super) fn value_error(error: impl std::fmt::Display) -> PyErr {
    crate::errors::ProviderError::new_err(error.to_string())
}

pub(super) fn index_error(error: impl std::fmt::Display) -> PyErr {
    PyIndexError::new_err(error.to_string())
}
