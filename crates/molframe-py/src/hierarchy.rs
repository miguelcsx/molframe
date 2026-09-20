//! Allocation-free Python navigation over the facade hierarchy.

use crate::bindings::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyAny;

fn position(index: isize, len: usize) -> PyResult<usize> {
    let len_signed = isize::try_from(len)
        .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("collection is too large"))?;
    let normalized = if index < 0 { len_signed + index } else { index };
    if normalized < 0 || normalized >= len_signed {
        return Err(pyo3::exceptions::PyIndexError::new_err(
            "hierarchy index is out of range",
        ));
    }
    usize::try_from(normalized)
        .map_err(|_| pyo3::exceptions::PyIndexError::new_err("invalid hierarchy index"))
}

fn contiguous_range<I>(mut indices: I) -> (u32, usize)
where
    I: Iterator<Item = u32>,
{
    match indices.next() {
        Some(first) => (first, 1 + indices.count()),
        None => (0, 0),
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Atom", frozen, skip_from_py_object)]
pub(crate) struct PyAtom {
    parent: PyStructure,
    index: u32,
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

    #[getter]
    fn coordinate(&self) -> Option<[f32; 3]> {
        self.parent
            .inner
            .coordinates()
            .get(self.index as usize)
            .copied()
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Atoms", frozen, skip_from_py_object)]
pub(crate) struct PyAtoms {
    parent: PyStructure,
    first: u32,
    len: usize,
}

#[pymethods]
impl PyAtoms {
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

    fn range(parent: PyStructure, first: u32, len: usize) -> Self {
        Self { parent, first, len }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Residue", frozen, skip_from_py_object)]
pub(crate) struct PyResidue {
    parent: PyStructure,
    index: u32,
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
}

#[derive(Clone, Debug)]
#[pyclass(name = "Residues", frozen, skip_from_py_object)]
pub(crate) struct PyResidues {
    parent: PyStructure,
    first: u32,
    len: usize,
}

#[pymethods]
impl PyResidues {
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

    fn range(parent: PyStructure, first: u32, len: usize) -> Self {
        Self { parent, first, len }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Chain", frozen, skip_from_py_object)]
pub(crate) struct PyChain {
    parent: PyStructure,
    index: u32,
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
    parent: PyStructure,
    first: u32,
    len: usize,
}

#[pymethods]
impl PyChains {
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
                .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(label.to_owned()))?
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

    fn range(parent: PyStructure, first: u32, len: usize) -> Self {
        Self { parent, first, len }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Model", frozen, skip_from_py_object)]
pub(crate) struct PyModel {
    parent: PyStructure,
    index: u32,
}

#[pymethods]
impl PyModel {
    #[getter]
    const fn index(&self) -> u32 {
        self.index
    }

    #[getter]
    fn number(&self) -> Option<i32> {
        self.parent
            .inner
            .models()
            .get(self.index as usize)
            .and_then(molframe::ModelRef::number)
    }

    #[getter]
    fn chains(&self) -> PyChains {
        let range = self
            .parent
            .inner
            .models()
            .get(self.index as usize)
            .map_or((0, 0), |model| {
                contiguous_range(model.chains().map(|chain| chain.index().get()))
            });
        PyChains::range(self.parent.clone(), range.0, range.1)
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Models", frozen, skip_from_py_object)]
pub(crate) struct PyModels {
    parent: PyStructure,
}

#[pymethods]
impl PyModels {
    fn __len__(&self) -> usize {
        self.parent.inner.model_count()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyModel> {
        let ordinal = position(index, self.parent.inner.model_count())?;
        let index = u32::try_from(ordinal)
            .map_err(|_| pyo3::exceptions::PyOverflowError::new_err("model index overflow"))?;
        Ok(PyModel {
            parent: self.parent.clone(),
            index,
        })
    }
}

impl PyModels {
    pub(crate) const fn new(parent: PyStructure) -> Self {
        Self { parent }
    }
}
