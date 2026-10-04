//! Models: one handle and the collection.

use super::{chain::PyChains, contiguous_range, position};
use crate::bindings::PyStructure;
use pyo3::prelude::*;

#[derive(Clone, Debug)]
#[pyclass(name = "Model", frozen, skip_from_py_object)]
pub(crate) struct PyModel {
    pub(super) parent: PyStructure,
    pub(super) index: u32,
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
    pub(super) parent: PyStructure,
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
