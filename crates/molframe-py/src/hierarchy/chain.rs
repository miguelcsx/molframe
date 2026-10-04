//! Chains: one handle and the views over a range.

use super::{contiguous_range, position, residue::PyResidues};
use crate::bindings::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyAny;

#[derive(Clone, Debug)]
#[pyclass(name = "Chain", frozen, skip_from_py_object)]
pub(crate) struct PyChain {
    pub(super) parent: PyStructure,
    pub(super) index: u32,
}

impl PyChain {
    fn handle(&self) -> Option<molframe::ChainRef<'_>> {
        self.parent.inner.chains().get(self.index as usize)
    }
}

#[pymethods]
impl PyChain {
    #[getter]
    const fn index(&self) -> u32 {
        self.index
    }

    #[getter]
    fn label(&self) -> Option<String> {
        self.parent
            .inner
            .chains()
            .get(self.index as usize)
            .and_then(|chain| chain.label().map(str::to_owned))
    }

    /// The depositor's chain label.
    #[getter]
    fn auth_label(&self) -> Option<String> {
        self.handle()
            .and_then(molframe::ChainRef::auth_label)
            .map(str::to_owned)
    }

    /// The index of the entity this chain instantiates.
    #[getter]
    fn entity(&self) -> Option<u32> {
        self.handle()
            .and_then(molframe::ChainRef::entity)
            .map(molframe::EntityIndex::get)
    }

    /// `"none"`, `"protein"`, `"dna"`, `"rna"`, `"nucleic_hybrid"`,
    /// `"saccharide"` or `"other"`.
    #[getter]
    fn polymer_kind(&self) -> Option<String> {
        self.handle()
            .map(|chain| crate::policy::snake(chain.polymer_kind().name()))
    }

    #[getter]
    fn residues(&self) -> PyResidues {
        let range = self
            .parent
            .inner
            .chains()
            .get(self.index as usize)
            .map_or((0, 0), |chain| {
                contiguous_range(chain.residues().map(|residue| residue.index().get()))
            });
        PyResidues::range(self.parent.clone(), range.0, range.1)
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Chains", frozen, skip_from_py_object)]
pub(crate) struct PyChains {
    pub(super) parent: PyStructure,
    pub(super) first: u32,
    pub(super) len: usize,
}

#[pymethods]
impl PyChains {
    /// The chains of this view as an Arrow stream, converted a batch at a time.
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, pyo3::types::PyCapsule>> {
        let rows = self.first as usize..self.first as usize + self.len;
        let table = molframe::interop::ChainTable::new(self.parent.inner.engine());
        crate::interop::arrow_capsule(py, table.arrow_stream_rows(rows))
    }
    fn __len__(&self) -> usize {
        self.len
    }

    fn __getitem__(&self, key: &Bound<'_, PyAny>) -> PyResult<PyChain> {
        let index = if let Ok(label) = key.extract::<&str>() {
            self.parent
                .inner
                .chains()
                .iter()
                .skip(self.first as usize)
                .take(self.len)
                .find(|chain| chain.label() == Some(label) || chain.auth_label() == Some(label))
                .map(|chain| chain.index().get())
                .ok_or_else(|| crate::error::key(label))?
        } else {
            let ordinal = position(key.extract::<isize>()?, self.len)?;
            let offset = u32::try_from(ordinal)
                .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("chain index overflow"))?;
            self.first
                .checked_add(offset)
                .ok_or_else(|| pyo3::exceptions::PyOverflowError::new_err("chain index overflow"))?
        };
        Ok(PyChain {
            parent: self.parent.clone(),
            index,
        })
    }
}

impl PyChains {
    pub(crate) fn all(parent: PyStructure) -> Self {
        let len = parent.inner.chain_count();
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
