//! Mechanical Python handles over the Rust hierarchy API.

use crate::atom::PyAtom;
use crate::errors::{index_error, key_error};
use crate::index::normalise_index;
use crate::structure::PyStructure;
use pdbiox::{AtomRef, ChainIndex, ModelIndex, ResidueIndex, Structure};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

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

#[pyclass(name = "Chains", frozen, skip_from_py_object)]
pub(crate) struct PyChains {
    inner: Structure,
    model: Option<ModelIndex>,
}

#[pyclass(name = "Chain", frozen, skip_from_py_object)]
pub(crate) struct PyChain {
    inner: Structure,
    index: ChainIndex,
}

#[pyclass(name = "Residues", frozen, skip_from_py_object)]
pub(crate) struct PyResidues {
    inner: Structure,
    chain: Option<ChainIndex>,
}

#[pyclass(name = "Residue", frozen, skip_from_py_object)]
pub(crate) struct PyResidue {
    inner: Structure,
    index: ResidueIndex,
}

#[pyclass(name = "ResidueAtoms", frozen, skip_from_py_object)]
pub(crate) struct PyResidueAtoms {
    inner: Structure,
    residue: ResidueIndex,
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
            .and_then(pdbiox::ModelRef::number)
    }

    #[getter]
    fn chains(&self) -> PyChains {
        PyChains {
            inner: self.inner.clone(),
            model: Some(self.local_index),
        }
    }
}

#[pymethods]
impl PyChains {
    fn __len__(&self) -> usize {
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

    fn chain_named(&self, py: Python<'_>, label: &str) -> PyResult<PyChain> {
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
    fn label(&self) -> Option<String> {
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
}

#[pymethods]
impl PyResidues {
    fn __len__(&self) -> usize {
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

#[pymethods]
impl PyResidue {
    #[getter]
    fn index(&self) -> u32 {
        self.index.get()
    }

    #[getter]
    fn name(&self) -> Option<String> {
        self.inner
            .residue(self.index)
            .and_then(|residue| residue.name().map(str::to_owned))
    }

    #[getter]
    fn label_seq_id(&self) -> Option<i32> {
        self.inner
            .residue(self.index)
            .and_then(pdbiox::ResidueRef::label_seq_id)
    }

    #[getter]
    fn auth_seq_id(&self) -> Option<i32> {
        self.inner
            .residue(self.index)
            .and_then(pdbiox::ResidueRef::auth_seq_id)
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
            .is_some_and(pdbiox::ResidueRef::is_het)
    }

    #[getter]
    fn atoms(&self) -> PyResidueAtoms {
        PyResidueAtoms {
            inner: self.inner.clone(),
            residue: self.index,
        }
    }
}

#[pymethods]
impl PyResidueAtoms {
    fn __len__(&self) -> usize {
        self.inner
            .residue(self.residue)
            .map_or(0, |residue| residue.atoms().count())
    }

    fn __getitem__(&self, key: &Bound<'_, PyAny>) -> PyResult<PyAtom> {
        let residue_key =
            isize::try_from(self.residue.get()).map_err(|_| index_error(key.py(), isize::MAX))?;
        let residue = self
            .inner
            .residue(self.residue)
            .ok_or_else(|| index_error(key.py(), residue_key))?;
        let atom = if let Ok(index) = key.extract::<isize>() {
            let Some(position) = normalise_index(index, residue.atoms().count()) else {
                return Err(index_error(key.py(), index));
            };
            residue.atom_at(position)
        } else if let Ok(name) = key.extract::<String>() {
            residue.atom(&name)
        } else {
            return Err(PyTypeError::new_err(
                "atom key must be an integer or string",
            ));
        };
        atom.map_or_else(
            || Err(key_error(key.py(), &key.str()?.to_string())),
            |atom| Ok(PyAtom::new(self.inner.clone(), AtomRef::index(atom))),
        )
    }
}

#[cfg(test)]
#[path = "hierarchy_tests.rs"]
mod tests;
