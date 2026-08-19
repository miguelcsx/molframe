//! Owned projections of structure data and coordinate-store variants.

use crate::bonds::PyBondTable;
use crate::core_annotations::PyAtomAnnotations;
use crate::core_records::PyAtomChunk;
use crate::core_storage::PyCoordinateBlock;
use crate::core_topology_root::PyTopology;
use crate::core_values::{PyCoordinateGeneration, PyInterner};
use crate::crystallography::PyUnitCell;
use crate::hierarchy::{
    PyChain, PyModel, PyResidue, atom_handle, chain_handle, model_handle, residue_handle,
};
use crate::index::PyModelIndex;
use crate::index::{PyAtomIndex, PyChainIndex, PyResidueIndex};
use crate::metadata::PyEntryMetadata;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "CoordinateStore", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCoordinateStore(pub(crate) pdbiox::CoordinateStore);

#[pyclass(name = "ExtensionStore", from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyExtensionStore(pub(crate) pdbiox::ExtensionStore);

#[pyclass(name = "StructureData", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructureData(pub(crate) pdbiox::StructureData);

#[pymethods]
impl PyCoordinateStore {
    #[staticmethod]
    fn single(block: PyCoordinateBlock) -> Self {
        Self(pdbiox::CoordinateStore::Single(block.0))
    }

    #[staticmethod]
    fn dense(frames: Vec<PyCoordinateBlock>) -> Self {
        Self(pdbiox::CoordinateStore::Dense {
            frames: frames.into_iter().map(|frame| frame.0).collect(),
        })
    }

    #[staticmethod]
    fn ragged(py: Python<'_>, models: Vec<Py<PyStructure>>) -> Self {
        Self(pdbiox::CoordinateStore::Ragged {
            models: models
                .into_iter()
                .map(|model| model.borrow(py).structure().clone())
                .collect(),
        })
    }

    fn model_count(&self) -> usize {
        self.0.model_count()
    }

    fn is_dense(&self) -> bool {
        self.0.is_dense()
    }

    fn block(&self, model: PyModelIndex) -> Option<PyCoordinateBlock> {
        self.0.block(model.0).cloned().map(PyCoordinateBlock)
    }

    fn ragged_models(&self) -> Option<Vec<PyStructure>> {
        self.0
            .ragged_models()
            .map(|models| models.iter().cloned().map(PyStructure::new).collect())
    }
}

#[pymethods]
impl PyExtensionStore {
    #[new]
    fn new() -> Self {
        Self(pdbiox::ExtensionStore::default())
    }

    fn clear(&mut self) {
        self.0.clear();
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn keys(&self) -> Vec<String> {
        self.0.keys().map(str::to_owned).collect()
    }
}

#[pymethods]
impl PyStructureData {
    #[staticmethod]
    fn empty() -> Self {
        Self(pdbiox::StructureData::empty())
    }

    fn atom_count(&self) -> u32 {
        self.0.chunks.last().map_or(0, |chunk| chunk.atoms().end)
    }

    fn model_count(&self) -> usize {
        self.0.coords.model_count()
    }

    #[getter]
    fn entry(&self) -> PyEntryMetadata {
        self.0.entry.clone().into()
    }

    #[getter]
    fn topology(&self) -> PyTopology {
        PyTopology(self.0.topology.clone())
    }

    #[getter]
    fn chunks(&self) -> Vec<PyAtomChunk> {
        self.0.chunks.iter().cloned().map(PyAtomChunk).collect()
    }

    #[getter]
    fn bonds(&self) -> PyBondTable {
        PyBondTable(self.0.bonds.clone())
    }

    #[getter]
    fn annotations(&self) -> PyAtomAnnotations {
        PyAtomAnnotations(self.0.annotations.clone())
    }

    #[getter]
    fn extensions(&self) -> PyExtensionStore {
        PyExtensionStore(self.0.extensions.clone())
    }

    #[getter]
    fn coords(&self) -> PyCoordinateStore {
        PyCoordinateStore(self.0.coords.clone())
    }

    #[getter]
    fn dictionary(&self) -> PyInterner {
        PyInterner(self.0.dictionary.clone())
    }

    #[getter]
    fn cell(&self) -> PyResult<Option<PyUnitCell>> {
        self.0.cell.map(PyUnitCell::from_native).transpose()
    }

    #[getter]
    fn generation(&self) -> PyCoordinateGeneration {
        PyCoordinateGeneration(self.0.generation)
    }

    fn to_structure(&self) -> PyStructure {
        PyStructure::new(pdbiox::Structure::new(self.0.clone()))
    }

    fn atom(&self, index: PyAtomIndex) -> Option<PyAtom> {
        let structure = self.to_structure();
        atom_handle(structure.structure(), index.0)
    }

    fn atoms(&self) -> Vec<PyAtom> {
        let structure = self.to_structure();
        (0..structure.structure().atom_count())
            .filter_map(|index| atom_handle(structure.structure(), pdbiox::AtomIndex::new(index)))
            .collect()
    }

    fn model(&self, index: PyModelIndex) -> Option<PyModel> {
        let structure = self.to_structure();
        model_handle(structure.structure(), index.0)
    }

    fn models(&self) -> Vec<PyModel> {
        let structure = self.to_structure();
        (0..structure.structure().model_count())
            .filter_map(|index| {
                u32::try_from(index).ok().and_then(|index| {
                    model_handle(structure.structure(), pdbiox::ModelIndex::new(index))
                })
            })
            .collect()
    }

    fn chain(&self, index: PyChainIndex) -> Option<PyChain> {
        let structure = self.to_structure();
        chain_handle(structure.structure(), index.0)
    }

    fn chains(&self) -> Vec<PyChain> {
        let structure = self.to_structure();
        (0..structure.structure().chain_count())
            .filter_map(|index| {
                u32::try_from(index).ok().and_then(|index| {
                    chain_handle(structure.structure(), pdbiox::ChainIndex::new(index))
                })
            })
            .collect()
    }

    fn chain_named(&self, label: &str) -> Option<PyChain> {
        let structure = self.to_structure();
        structure
            .structure()
            .data()
            .chain_named(label)
            .and_then(|chain| chain_handle(structure.structure(), chain.index()))
    }

    fn residue(&self, index: PyResidueIndex) -> Option<PyResidue> {
        let structure = self.to_structure();
        residue_handle(structure.structure(), index.0)
    }

    fn residues(&self) -> Vec<PyResidue> {
        let structure = self.to_structure();
        (0..structure.structure().residue_count())
            .filter_map(|index| {
                u32::try_from(index).ok().and_then(|index| {
                    residue_handle(structure.structure(), pdbiox::ResidueIndex::new(index))
                })
            })
            .collect()
    }
}

#[pymethods]
impl PyStructure {
    fn data(&self) -> PyStructureData {
        PyStructureData(self.structure().data().clone())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCoordinateStore>()?;
    module.add_class::<PyExtensionStore>()?;
    module.add_class::<PyStructureData>()?;
    Ok(())
}
use crate::atom::PyAtom;
