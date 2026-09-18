//! Mechanical Python handles over the Rust hierarchy API.
//!
//! Every handle is a structure snapshot plus an index, so a handle reached
//! through the full table and one reached through its parent address the same
//! rows. The three sibling modules extend the same classes — atoms, chains and
//! residues — rather than defining parallel ones.

use crate::atom::PyAtom;
use crate::errors::{index_error, key_error};
use crate::index::{PyAtomIndex, PyChainIndex, PyModelIndex, PyResidueIndex, normalise_index};
use crate::structure::PyStructure;
use molframe::{AtomRef, ChainIndex, ModelIndex, ResidueIndex, Structure};
use pyo3::prelude::*;

#[path = "hierarchy_atoms.rs"]
mod hierarchy_atoms;
#[path = "hierarchy_chains.rs"]
mod hierarchy_chains;
#[path = "hierarchy_residues.rs"]
mod hierarchy_residues;

pub(crate) use hierarchy_chains::{PyChain, PyChains};
pub(crate) use hierarchy_residues::{PyResidue, PyResidues};

#[pyclass(name = "Models", frozen, skip_from_py_object)]
pub(crate) struct PyModels {
    inner: Structure,
}

#[pyclass(name = "Model", frozen, skip_from_py_object)]
pub(crate) struct PyModel {
    inner: Structure,
    index: ModelIndex,
    local_index: ModelIndex,
}

#[pyclass(name = "ResidueAtoms", frozen, skip_from_py_object)]
pub(crate) struct PyResidueAtoms {
    inner: Structure,
    residue: ResidueIndex,
}

pub(crate) fn atom_handle(structure: &Structure, index: molframe::AtomIndex) -> Option<PyAtom> {
    structure
        .atom(index)
        .map(|atom| PyAtom::new(structure.clone(), AtomRef::index(atom)))
}

pub(crate) fn model_handle(structure: &Structure, index: ModelIndex) -> Option<PyModel> {
    structure
        .model_snapshot(index)
        .map(|(inner, local_index)| PyModel {
            inner,
            index,
            local_index,
        })
}

pub(crate) fn chain_handle(structure: &Structure, index: ChainIndex) -> Option<PyChain> {
    structure.chain(index).map(|_| PyChain {
        inner: structure.clone(),
        index,
    })
}

pub(crate) fn residue_handle(structure: &Structure, index: ResidueIndex) -> Option<PyResidue> {
    structure.residue(index).map(|_| PyResidue {
        inner: structure.clone(),
        index,
    })
}

#[pymethods]
impl PyStructure {
    #[getter]
    fn models(&self) -> PyModels {
        PyModels {
            inner: self.structure().clone(),
        }
    }

    #[getter]
    fn chains(&self) -> PyChains {
        PyChains::new_full(self.structure().clone())
    }

    #[getter]
    fn residues(&self) -> PyResidues {
        PyResidues::new_full(self.structure().clone())
    }

    fn atom(&self, index: PyAtomIndex) -> Option<PyAtom> {
        atom_handle(self.structure(), index.0)
    }

    fn model(&self, index: PyModelIndex) -> Option<PyModel> {
        model_handle(self.structure(), index.0)
    }

    fn chain(&self, index: PyChainIndex) -> Option<PyChain> {
        chain_handle(self.structure(), index.0)
    }

    fn residue(&self, index: PyResidueIndex) -> Option<PyResidue> {
        residue_handle(self.structure(), index.0)
    }
}

#[pymethods]
impl PyModels {
    fn __len__(&self) -> usize {
        self.inner.model_count()
    }

    fn __getitem__(&self, py: Python<'_>, python_index: isize) -> PyResult<PyModel> {
        let Some(position) = normalise_index(python_index, self.inner.model_count()) else {
            return Err(index_error(py, python_index));
        };
        let position = u32::try_from(position).map_err(|_| index_error(py, python_index))?;
        let index = ModelIndex::new(position);
        let Some((inner, local_index)) = self.inner.model_snapshot(index) else {
            return Err(index_error(py, python_index));
        };
        Ok(PyModel {
            inner,
            index,
            local_index,
        })
    }
}

#[pymethods]
impl PyModel {
    #[getter]
    fn index(&self) -> u32 {
        self.index.get()
    }

    #[getter]
    fn number(&self) -> Option<i32> {
        self.inner
            .model(self.local_index)
            .and_then(molframe::ModelRef::number)
    }

    #[getter]
    fn index_id(&self) -> PyModelIndex {
        PyModelIndex(self.index)
    }

    fn chain_at(&self, py: Python<'_>, position: isize) -> PyResult<PyChain> {
        let Some(position) = normalise_index(position, self.chains().__len__()) else {
            return Err(index_error(py, position));
        };
        self.inner
            .model(self.local_index)
            .and_then(|model| model.chain_at(position))
            .map(|chain| PyChain {
                inner: self.inner.clone(),
                index: chain.index(),
            })
            .ok_or_else(|| index_error(py, position.cast_signed()))
    }

    fn chain(&self, py: Python<'_>, label: &str) -> PyResult<PyChain> {
        self.inner
            .model(self.local_index)
            .and_then(|model| model.chain(label))
            .map(|chain| PyChain {
                inner: self.inner.clone(),
                index: chain.index(),
            })
            .ok_or_else(|| key_error(py, label))
    }

    #[getter]
    fn chains(&self) -> PyChains {
        PyChains {
            inner: self.inner.clone(),
            model: Some(self.local_index),
        }
    }
}

#[cfg(test)]
#[path = "hierarchy_tests.rs"]
mod tests;
