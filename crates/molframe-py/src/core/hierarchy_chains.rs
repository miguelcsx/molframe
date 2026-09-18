//! Chain-scoped Python handles.
//!
//! A [`PyChains`] is either the whole table or one model's chains; both answer
//! `len` and both resolve an integer or a label to a [`PyChain`]. The handles
//! hold a structure snapshot and an index, so a chain reached through a model
//! and one reached through the full table address the same rows.

use crate::core_values::PySymbolId;
use crate::errors::{index_error, key_error};
use crate::index::{PyChainIndex, normalise_index};
use molframe::{ChainIndex, ModelIndex, Structure};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

use super::hierarchy_residues::{PyResidue, PyResidues};

#[pyclass(name = "Chains", frozen, skip_from_py_object)]
pub(crate) struct PyChains {
    pub(super) inner: Structure,
    pub(super) model: Option<ModelIndex>,
}

#[pyclass(name = "Chain", frozen, skip_from_py_object)]
pub(crate) struct PyChain {
    pub(crate) inner: Structure,
    pub(crate) index: ChainIndex,
}

#[pymethods]
impl PyChains {
    pub(super) fn __len__(&self) -> usize {
        match self.model.and_then(|index| self.inner.model(index)) {
            Some(model) => model.chains().count(),
            None => self.inner.chain_count(),
        }
    }

    fn __getitem__(&self, key: &Bound<'_, PyAny>) -> PyResult<PyChain> {
        if let Ok(index) = key.extract::<isize>() {
            return self.chain_at(key.py(), index);
        }
        if let Ok(label) = key.extract::<String>() {
            return self.chain_named(key.py(), &label);
        }
        Err(PyTypeError::new_err(
            "chain key must be an integer or string",
        ))
    }
}

impl PyChains {
    pub(crate) const fn new_full(inner: Structure) -> Self {
        Self { inner, model: None }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }

    pub(crate) const fn is_full_table(&self) -> bool {
        self.model.is_none()
    }

    fn chain_at(&self, py: Python<'_>, index: isize) -> PyResult<PyChain> {
        let Some(position) = normalise_index(index, self.__len__()) else {
            return Err(index_error(py, index));
        };
        let found = if let Some(model) = self.model.and_then(|model| self.inner.model(model)) {
            model.chain_at(position)
        } else {
            let position = u32::try_from(position).map_err(|_| index_error(py, index))?;
            self.inner.chain(ChainIndex::new(position))
        };
        found.map_or_else(
            || Err(index_error(py, index)),
            |chain| {
                Ok(PyChain {
                    inner: self.inner.clone(),
                    index: chain.index(),
                })
            },
        )
    }

    pub(super) fn chain_named(&self, py: Python<'_>, label: &str) -> PyResult<PyChain> {
        let found = match self.model.and_then(|model| self.inner.model(model)) {
            Some(model) => model.chain(label),
            None => self.inner.data().chain_named(label),
        };
        found.map_or_else(
            || Err(key_error(py, label)),
            |chain| {
                Ok(PyChain {
                    inner: self.inner.clone(),
                    index: chain.index(),
                })
            },
        )
    }
}

#[pymethods]
impl PyChain {
    #[getter]
    fn index(&self) -> u32 {
        self.index.get()
    }

    #[getter]
    fn index_id(&self) -> PyChainIndex {
        PyChainIndex(self.index)
    }

    #[getter]
    fn label_asym_id(&self) -> Option<PySymbolId> {
        self.inner
            .chain(self.index)
            .and_then(molframe::ChainRef::label_asym_id)
            .map(PySymbolId)
    }

    #[getter]
    fn auth_asym_id(&self) -> Option<PySymbolId> {
        self.inner
            .chain(self.index)
            .and_then(molframe::ChainRef::auth_asym_id)
            .map(PySymbolId)
    }

    #[getter]
    pub(super) fn label(&self) -> Option<String> {
        self.inner
            .chain(self.index)
            .and_then(|chain| chain.label().map(str::to_owned))
    }

    #[getter]
    fn auth_label(&self) -> Option<String> {
        self.inner
            .chain(self.index)
            .and_then(|chain| chain.auth_label().map(str::to_owned))
    }

    #[getter]
    fn residues(&self) -> PyResidues {
        PyResidues {
            inner: self.inner.clone(),
            chain: Some(self.index),
        }
    }

    fn residue_at(&self, py: Python<'_>, position: isize) -> PyResult<PyResidue> {
        let Some(position) = normalise_index(position, self.residues().__len__()) else {
            return Err(index_error(py, position));
        };
        self.inner
            .chain(self.index)
            .and_then(|chain| chain.residue_at(position))
            .map(|residue| PyResidue {
                inner: self.inner.clone(),
                index: residue.index(),
            })
            .ok_or_else(|| index_error(py, position.cast_signed()))
    }

    fn residue(&self, number: i32) -> Option<PyResidue> {
        self.inner
            .chain(self.index)
            .and_then(|chain| chain.residue(number))
            .map(|residue| PyResidue {
                inner: self.inner.clone(),
                index: residue.index(),
            })
    }
}
