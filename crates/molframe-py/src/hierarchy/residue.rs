//! Residues: one handle, a selection of them, and the views over a range.

use super::{
    atom::{PyAtom, PyAtoms},
    contiguous_range, position,
};
use crate::bindings::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyAny;

#[derive(Clone, Debug)]
#[pyclass(name = "Residue", frozen, skip_from_py_object)]
pub(crate) struct PyResidue {
    pub(super) parent: PyStructure,
    pub(super) index: u32,
}

impl PyResidue {
    fn handle(&self) -> Option<molframe::ResidueRef<'_>> {
        self.parent.inner.residues().get(self.index as usize)
    }
}

#[pymethods]
impl PyResidue {
    #[getter]
    const fn index(&self) -> u32 {
        self.index
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.parent
            .inner
            .residues()
            .get(self.index as usize)
            .and_then(|residue| residue.name().map(str::to_owned))
    }

    /// The sequence number in the file's normalised numbering.
    #[getter]
    fn number(&self) -> Option<i32> {
        self.handle().and_then(molframe::ResidueRef::label_seq_id)
    }

    /// The depositor's sequence number.
    #[getter]
    fn auth_number(&self) -> Option<i32> {
        self.handle().and_then(molframe::ResidueRef::auth_seq_id)
    }

    /// The depositor's residue name.
    #[getter]
    fn auth_name(&self) -> Option<String> {
        self.handle()
            .and_then(molframe::ResidueRef::auth_name)
            .map(str::to_owned)
    }

    /// The insertion code, or `None` for a residue without one.
    #[getter]
    fn insertion_code(&self) -> Option<String> {
        self.handle()
            .and_then(molframe::ResidueRef::ins_code)
            .map(str::to_owned)
    }

    /// Whether the file lists the residue as a heterogen.
    #[getter]
    fn is_hetero(&self) -> bool {
        self.handle().is_some_and(molframe::ResidueRef::is_het)
    }

    /// The chain the residue belongs to.
    #[getter]
    fn chain(&self) -> Option<super::chain::PyChain> {
        self.parent
            .inner
            .engine()
            .data()
            .topology
            .chains
            .containing(self.index)
            .map(|chain| super::chain::PyChain {
                parent: self.parent.clone(),
                index: chain.get(),
            })
    }

    /// The residue's secondary-structure state.
    #[getter]
    fn secondary_structure(&self) -> Option<crate::secondary::PySecondaryStructure> {
        self.parent
            .inner
            .secondary_structure()
            .get(self.index as usize)
            .copied()
            .map(Into::into)
    }

    /// Which assigner produced the residue's secondary-structure state.
    #[getter]
    fn secondary_source(&self) -> Option<crate::secondary::PySecondarySource> {
        self.parent
            .inner
            .secondary_source()
            .get(self.index as usize)
            .copied()
            .map(Into::into)
    }

    #[getter]
    fn atoms(&self) -> PyAtoms {
        let range = self
            .parent
            .inner
            .residues()
            .get(self.index as usize)
            .map_or((0, 0), |residue| {
                contiguous_range(residue.atoms().map(|atom| atom.index().get()))
            });
        PyAtoms::range(self.parent.clone(), range.0, range.1)
    }
    fn atom(&self, name: &str) -> Option<PyAtom> {
        self.parent
            .inner
            .residues()
            .get(self.index as usize)
            .and_then(|residue| residue.atom(name))
            .map(|atom| PyAtom {
                parent: self.parent.clone(),
                index: atom.index().get(),
            })
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "ResidueSelection", frozen, skip_from_py_object)]
pub(crate) struct PyResidueSelection {
    pub(super) parent: PyStructure,
    pub(super) indices: Vec<u32>,
}

#[pymethods]
impl PyResidueSelection {
    const fn __len__(&self) -> usize {
        self.indices.len()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyResidue> {
        let index = position(index, self.indices.len())?;
        Ok(PyResidue {
            parent: self.parent.clone(),
            index: self.indices[index],
        })
    }
}

impl PyResidueSelection {
    pub(crate) fn from_selection(parent: PyStructure, selection: &molframe::Selection) -> Self {
        let mut indices = selection
            .atoms()
            .filter_map(|atom| atom.residue().map(|residue| residue.index().get()))
            .collect::<Vec<_>>();
        indices.sort_unstable();
        indices.dedup();
        Self { parent, indices }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Residues", frozen, skip_from_py_object)]
pub(crate) struct PyResidues {
    pub(super) parent: PyStructure,
    pub(super) first: u32,
    pub(super) len: usize,
}

#[pymethods]
impl PyResidues {
    /// The residues of this view as an Arrow stream, converted a batch at a time.
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, pyo3::types::PyCapsule>> {
        let rows = self.first as usize..self.first as usize + self.len;
        let table = molframe::interop::ResidueTable::new(self.parent.inner.engine());
        crate::interop::arrow_capsule(py, table.arrow_stream_rows(rows))
    }
    const fn __len__(&self) -> usize {
        self.len
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyResidue> {
        let offset = position(index, self.len)?;
        let offset = u32::try_from(offset)
            .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("residue index overflow"))?;
        let index = self
            .first
            .checked_add(offset)
            .ok_or_else(|| pyo3::exceptions::PyOverflowError::new_err("residue index overflow"))?;
        Ok(PyResidue {
            parent: self.parent.clone(),
            index,
        })
    }
}

impl PyResidues {
    pub(crate) fn all(parent: PyStructure) -> Self {
        let len = parent.inner.residue_count();
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
