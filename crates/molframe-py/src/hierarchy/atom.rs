//! Atoms: one handle, and the views over a range of them.

use super::{position, residue::PyResidue};
use crate::bindings::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyAny;

#[derive(Clone, Debug)]
#[pyclass(name = "Atom", frozen, skip_from_py_object)]
pub(crate) struct PyAtom {
    pub(super) parent: PyStructure,
    pub(super) index: u32,
}

impl PyAtom {
    fn handle(&self) -> Option<molframe::AtomRef<'_>> {
        self.parent.inner.atoms().get(self.index as usize)
    }
}

#[pymethods]
impl PyAtom {
    #[getter]
    const fn index(&self) -> u32 {
        self.index
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.parent
            .inner
            .atoms()
            .get(self.index as usize)
            .and_then(|atom| atom.name().map(str::to_owned))
    }

    /// The effective chemical component identifier for this atom.
    ///
    /// This is the native component identity, including an atom-level
    /// alternate component when one is present; no chemistry is inferred.
    #[getter]
    fn component_name(&self) -> Option<String> {
        self.parent
            .inner
            .atoms()
            .get(self.index as usize)
            .and_then(molframe::AtomRef::component_name)
            .map(str::to_owned)
    }
    /// The element symbol, or `None` for an unknown element.
    #[getter]
    fn element(&self) -> Option<String> {
        self.handle()
            .and_then(molframe::AtomRef::element)
            .map(|element| element.symbol().to_owned())
    }

    /// The atomic number, or `None` for an unknown element.
    #[getter]
    fn atomic_number(&self) -> Option<u8> {
        self.handle()
            .and_then(molframe::AtomRef::element)
            .map(molframe::Element::atomic_number)
    }

    /// Occupancy, or `None` where the file records none.
    #[getter]
    fn occupancy(&self) -> Option<f32> {
        self.handle().and_then(molframe::AtomRef::occupancy)
    }

    /// Temperature factor, or `None` where the file records none.
    #[getter]
    fn b_factor(&self) -> Option<f32> {
        self.handle().and_then(molframe::AtomRef::b_factor)
    }

    /// Formal charge, or `None` where none is recorded.
    #[getter]
    fn formal_charge(&self) -> Option<i8> {
        self.handle().and_then(molframe::AtomRef::formal_charge)
    }

    /// The alternate-location label, or `None` for an atom in every conformer.
    #[getter]
    fn altloc(&self) -> Option<String> {
        self.handle()
            .and_then(molframe::AtomRef::alt_label)
            .map(str::to_owned)
    }

    /// The serial number the file gave the atom.
    #[getter]
    fn serial(&self) -> Option<u32> {
        self.handle().and_then(molframe::AtomRef::atom_site_id)
    }

    /// The depositor's atom name, where it differs from the normalised one.
    #[getter]
    fn auth_name(&self) -> Option<String> {
        self.handle()
            .and_then(molframe::AtomRef::auth_name)
            .map(str::to_owned)
    }

    /// The six anisotropic displacement parameters `[U11, U22, U33, U12, U13, U23]`
    /// in Å², or `None` for an isotropic atom.
    #[getter]
    fn anisotropy(&self) -> Option<[f32; 6]> {
        self.parent
            .inner
            .anisotropy()
            .for_atom(molframe::AtomIndex::new(self.index))
    }

    #[getter]
    fn coordinate(&self) -> Option<[f32; 3]> {
        self.parent
            .inner
            .coordinates()
            .get(self.index as usize)
            .copied()
    }
    #[getter]
    fn residue(&self) -> Option<PyResidue> {
        self.parent
            .inner
            .atoms()
            .get(self.index as usize)
            .and_then(molframe::AtomRef::residue)
            .map(|residue| PyResidue {
                parent: self.parent.clone(),
                index: residue.index().get(),
            })
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Atoms", frozen, skip_from_py_object)]
pub(crate) struct PyAtoms {
    pub(super) parent: PyStructure,
    pub(super) first: u32,
    pub(super) len: usize,
}

#[pymethods]
impl PyAtoms {
    /// The atoms of this view as an Arrow stream, converted a batch at a time.
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, pyo3::types::PyCapsule>> {
        let rows = self.first as usize..self.first as usize + self.len;
        let table = molframe::interop::AtomTable::new(self.parent.inner.engine());
        crate::interop::arrow_capsule(py, table.arrow_stream_rows(rows))
    }
    const fn __len__(&self) -> usize {
        self.len
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyAtom> {
        let offset = position(index, self.len)?;
        let offset = u32::try_from(offset)
            .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("atom index overflow"))?;
        let index = self
            .first
            .checked_add(offset)
            .ok_or_else(|| pyo3::exceptions::PyOverflowError::new_err("atom index overflow"))?;
        Ok(PyAtom {
            parent: self.parent.clone(),
            index,
        })
    }
}

impl PyAtoms {
    pub(crate) fn all(parent: PyStructure) -> Self {
        let len = parent.inner.atom_count() as usize;
        Self {
            parent,
            first: 0,
            len,
        }
    }

    pub(super) fn range(parent: PyStructure, first: u32, len: usize) -> Self {
        Self { parent, first, len }
    }
}
