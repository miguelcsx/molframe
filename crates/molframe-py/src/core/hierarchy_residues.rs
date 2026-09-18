//! Residue-scoped Python handles.
//!
//! A [`PyResidues`] is either the whole table or one chain's residues, and a
//! [`PyResidue`] exposes both identifier namespaces the format may carry. The
//! handles hold a structure snapshot and an index, so nothing here copies a
//! coordinate or re-reads the source file.

use crate::atom::PyAtom;
use crate::core_values::PySymbolId;
use crate::errors::{index_error, key_error};
use crate::index::{PyResidueIndex, normalise_index};
use molframe::{AtomRef, ChainIndex, ResidueIndex, Structure};
use pyo3::prelude::*;

use super::PyResidueAtoms;

#[pyclass(name = "Residues", frozen, skip_from_py_object)]
pub(crate) struct PyResidues {
    pub(super) inner: Structure,
    pub(super) chain: Option<ChainIndex>,
}

#[pyclass(name = "Residue", frozen, skip_from_py_object)]
pub(crate) struct PyResidue {
    pub(super) inner: Structure,
    pub(super) index: ResidueIndex,
}

#[pymethods]
impl PyResidues {
    pub(super) fn __len__(&self) -> usize {
        match self.chain.and_then(|index| self.inner.chain(index)) {
            Some(chain) => chain.residues().count(),
            None => self.inner.residue_count(),
        }
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<PyResidue> {
        let Some(position) = normalise_index(index, self.__len__()) else {
            return Err(index_error(py, index));
        };
        let found = if let Some(chain) = self.chain.and_then(|chain| self.inner.chain(chain)) {
            chain.residue_at(position)
        } else {
            let position = u32::try_from(position).map_err(|_| index_error(py, index))?;
            self.inner.residue(ResidueIndex::new(position))
        };
        found.map_or_else(
            || Err(index_error(py, index)),
            |residue| {
                Ok(PyResidue {
                    inner: self.inner.clone(),
                    index: residue.index(),
                })
            },
        )
    }
}

impl PyResidues {
    pub(crate) const fn new_full(inner: Structure) -> Self {
        Self { inner, chain: None }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }

    pub(crate) const fn is_full_table(&self) -> bool {
        self.chain.is_none()
    }
}

impl PyResidue {
    pub(crate) const fn new(inner: Structure, index: ResidueIndex) -> Self {
        Self { inner, index }
    }
}

#[pymethods]
impl PyResidue {
    #[getter]
    fn index(&self) -> u32 {
        self.index.get()
    }

    #[getter]
    fn index_id(&self) -> PyResidueIndex {
        PyResidueIndex(self.index)
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner
            .residue(self.index)
            .and_then(|residue| residue.name().map(str::to_owned))
    }

    #[getter]
    fn label_comp_id(&self) -> Option<PySymbolId> {
        self.inner
            .residue(self.index)
            .and_then(molframe::ResidueRef::label_comp_id)
            .map(PySymbolId)
    }

    #[getter]
    fn auth_comp_id(&self) -> Option<PySymbolId> {
        self.inner
            .residue(self.index)
            .and_then(molframe::ResidueRef::auth_comp_id)
            .map(PySymbolId)
    }

    #[getter]
    fn auth_name(&self) -> Option<String> {
        self.inner
            .residue(self.index)
            .and_then(|residue| residue.auth_name().map(str::to_owned))
    }

    #[getter]
    fn label_seq_id(&self) -> Option<i32> {
        self.inner
            .residue(self.index)
            .and_then(molframe::ResidueRef::label_seq_id)
    }

    #[getter]
    fn auth_seq_id(&self) -> Option<i32> {
        self.inner
            .residue(self.index)
            .and_then(molframe::ResidueRef::auth_seq_id)
    }

    #[getter]
    fn ins_code(&self) -> Option<String> {
        self.inner
            .residue(self.index)
            .and_then(|residue| residue.ins_code().map(str::to_owned))
    }

    #[getter]
    fn is_het(&self) -> bool {
        self.inner
            .residue(self.index)
            .is_some_and(molframe::ResidueRef::is_het)
    }

    #[getter]
    fn atoms(&self) -> PyResidueAtoms {
        PyResidueAtoms {
            inner: self.inner.clone(),
            residue: self.index,
        }
    }

    fn atom_at(&self, py: Python<'_>, position: isize) -> PyResult<PyAtom> {
        let atom_count = self
            .inner
            .residue(self.index)
            .map_or(0, |residue| residue.atoms().count());
        let Some(position) = normalise_index(position, atom_count) else {
            return Err(index_error(py, position));
        };
        self.inner
            .residue(self.index)
            .and_then(|residue| residue.atom_at(position))
            .map(|atom| PyAtom::new(self.inner.clone(), AtomRef::index(atom)))
            .ok_or_else(|| index_error(py, position.cast_signed()))
    }

    fn atom(&self, py: Python<'_>, name: &str) -> PyResult<PyAtom> {
        self.inner
            .residue(self.index)
            .and_then(|residue| residue.atom(name))
            .map(|atom| PyAtom::new(self.inner.clone(), AtomRef::index(atom)))
            .ok_or_else(|| key_error(py, name))
    }
}
