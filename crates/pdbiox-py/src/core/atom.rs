//! Lightweight Python atom handles over one shared Rust snapshot.

use crate::chemistry::PyElement;
use crate::core_values::{PyAltId, PySymbolId};
use crate::errors::index_error;
use crate::hierarchy::PyResidue;
use crate::index::{PyAtomIndex, normalise_index};
use crate::structure::PyStructure;
use pdbiox::{AtomIndex, AtomRef, Structure};
use pyo3::prelude::*;

#[pyclass(name = "Atoms", frozen, skip_from_py_object)]
pub(crate) struct PyAtoms {
    inner: Structure,
}

#[pyclass(name = "Atom", frozen, skip_from_py_object)]
pub(crate) struct PyAtom {
    inner: Structure,
    index: AtomIndex,
}

impl PyAtoms {
    pub(crate) const fn new(inner: Structure) -> Self {
        Self { inner }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }
}

impl PyAtom {
    pub(crate) const fn new(inner: Structure, index: AtomIndex) -> Self {
        Self { inner, index }
    }
}

#[pymethods]
impl PyAtom {
    fn alt_label(&self) -> Option<String> {
        self.altloc()
    }

    fn component_name(&self) -> Option<String> {
        self.component()
    }

    fn position(&self) -> Option<[f32; 3]> {
        self.coord()
    }
}

#[pymethods]
impl PyStructure {
    #[getter]
    fn atoms(&self) -> PyAtoms {
        PyAtoms::new(self.structure().clone())
    }
}

#[pymethods]
impl PyAtoms {
    fn __len__(&self) -> usize {
        self.inner.atom_count() as usize
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<PyAtom> {
        let count = self.inner.atom_count() as usize;
        let Some(position) = normalise_index(index, count) else {
            return Err(index_error(py, index));
        };
        let position = u32::try_from(position).map_err(|_| index_error(py, index))?;
        Ok(PyAtom::new(self.inner.clone(), AtomIndex::new(position)))
    }
}

#[pymethods]
impl PyAtom {
    #[getter]
    fn index(&self) -> u32 {
        self.index.get()
    }

    #[getter]
    fn index_id(&self) -> PyAtomIndex {
        PyAtomIndex(self.index)
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::name)
            .map(str::to_owned)
    }

    #[getter]
    fn name_symbol(&self) -> Option<PySymbolId> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::name_symbol)
            .map(PySymbolId)
    }

    #[getter]
    fn auth_name(&self) -> Option<String> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::auth_name)
            .map(str::to_owned)
    }

    #[getter]
    fn auth_name_symbol(&self) -> Option<PySymbolId> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::auth_name_symbol)
            .map(PySymbolId)
    }

    #[getter]
    fn altloc(&self) -> Option<String> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::alt_label)
            .map(str::to_owned)
    }

    #[getter]
    fn alt_id(&self) -> Option<PyAltId> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::alt_id)
            .map(PyAltId)
    }

    #[getter]
    fn component(&self) -> Option<String> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::component_name)
            .map(str::to_owned)
    }

    #[getter]
    fn component_id(&self) -> Option<PySymbolId> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::component_id)
            .map(PySymbolId)
    }

    #[getter]
    fn element(&self) -> Option<&'static str> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::element)
            .map(pdbiox::Element::symbol)
    }

    #[getter]
    fn element_value(&self) -> Option<PyElement> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::element)
            .map(PyElement)
    }

    #[getter]
    fn coord(&self) -> Option<[f32; 3]> {
        self.inner.atom(self.index).and_then(AtomRef::position)
    }

    #[getter]
    fn b_factor(&self) -> Option<f32> {
        self.inner.atom(self.index).and_then(AtomRef::b_factor)
    }

    #[getter]
    fn occupancy(&self) -> Option<f32> {
        self.inner.atom(self.index).and_then(AtomRef::occupancy)
    }

    #[getter]
    fn formal_charge(&self) -> Option<i8> {
        self.inner.atom(self.index).and_then(AtomRef::formal_charge)
    }

    #[getter]
    fn atom_site_id(&self) -> Option<u32> {
        self.inner.atom(self.index).and_then(AtomRef::atom_site_id)
    }

    #[getter]
    fn residue(&self) -> Option<PyResidue> {
        self.inner.atom(self.index).and_then(|atom| {
            atom.residue()
                .map(|residue| PyResidue::new(self.inner.clone(), residue.index()))
        })
    }
}
