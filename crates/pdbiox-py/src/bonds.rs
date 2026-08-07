//! Mechanical Python handle over the Rust bond table.

use crate::structure::PyStructure;
use pdbiox::Structure;
use pyo3::prelude::*;

#[pyclass(name = "Bonds", frozen, skip_from_py_object)]
pub(crate) struct PyBonds {
    inner: Structure,
}

impl PyBonds {
    pub(crate) const fn new(inner: Structure) -> Self {
        Self { inner }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }
}

#[pymethods]
impl PyStructure {
    #[getter]
    fn bonds(&self) -> PyBonds {
        PyBonds::new(self.structure().clone())
    }
}

#[pymethods]
impl PyBonds {
    fn __len__(&self) -> usize {
        self.inner.data().bonds.len()
    }
}
