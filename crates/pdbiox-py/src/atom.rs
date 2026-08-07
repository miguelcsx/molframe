//! Lightweight Python atom handles over one shared Rust snapshot.

use crate::errors::index_error;
use crate::index::normalise_index;
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
    fn name(&self) -> Option<String> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::name)
            .map(str::to_owned)
    }

    #[getter]
    fn element(&self) -> Option<&'static str> {
        self.inner
            .atom(self.index)
            .and_then(AtomRef::element)
            .map(pdbiox::Element::symbol)
    }

    #[getter]
    fn coord(&self) -> Option<[f32; 3]> {
        self.inner.atom(self.index).and_then(AtomRef::position)
    }

    #[getter]
    fn b_factor(&self) -> Option<f32> {
        self.inner.atom(self.index).and_then(AtomRef::b_factor)
    }
}
