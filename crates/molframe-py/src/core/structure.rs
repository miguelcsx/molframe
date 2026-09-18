//! Python handles and zero-copy array conversion for a structure snapshot.

use crate::contract::PyDiagnostic;
use crate::core_annotations::PyAtomAnnotations;
use crate::core_data::PyExtensionStore;
use crate::core_values::{PyCoordinateGeneration, PySymbolId};
use crate::errors::read_error;
use crate::index::PyModelIndex;
use molframe::Structure;
use numpy::ndarray::ArrayView2;
use numpy::{PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use std::path::PathBuf;

#[pyclass(name = "Structure", skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructure {
    inner: Structure,
}

#[pyclass(frozen)]
struct SnapshotOwner {
    inner: Structure,
}

impl PyStructure {
    pub(crate) const fn new(inner: Structure) -> Self {
        Self { inner }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }

    pub(crate) fn replace(&mut self, inner: Structure) {
        self.inner = inner;
    }

    fn readonly_positions(
        py: Python<'_>,
        inner: Structure,
        model: Option<molframe::ModelIndex>,
    ) -> PyResult<Bound<'_, PyArray2<f32>>> {
        let owner = Bound::new(py, SnapshotOwner { inner })?;
        let (position_count, pointer) = {
            let snapshot = owner.borrow();
            let positions = match model.and_then(|index| snapshot.inner.model_positions(index)) {
                Some(positions) => positions,
                None => snapshot.inner.positions(),
            };
            (positions.len(), positions.as_ptr().cast::<f32>())
        };
        position_count.checked_mul(3).ok_or_else(|| {
            pyo3::exceptions::PyOverflowError::new_err("coordinate array shape overflows usize")
        })?;
        let view = unsafe { ArrayView2::from_shape_ptr((position_count, 3), pointer) };
        let array = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
        array.readwrite().make_nonwriteable();
        Ok(array)
    }
}

#[pymethods]
impl PyStructure {
    #[staticmethod]
    fn empty() -> Self {
        Self::new(molframe::Structure::new(molframe::StructureData::empty()))
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    /// Publishes the immutable native snapshot for sibling Rust extensions.
    ///
    /// The capsule is deliberately private to the Python protocol: consumers
    /// must validate its name before reading the `repr(C)` payload.
    fn _pdviewx_structure_capsule<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyCapsule>> {
        PyCapsule::new_with_value(py, self.inner.clone(), c"molframe.Structure")
    }

    fn write(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        molframe::write(path, &self.inner).map_err(|findings| read_error(py, &findings))
    }

    #[getter]
    fn atom_count(&self) -> u32 {
        self.inner.atom_count()
    }

    #[getter]
    fn model_count(&self) -> usize {
        self.inner.model_count()
    }

    #[getter]
    fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    #[getter]
    fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    #[getter]
    fn xyz<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        Self::readonly_positions(py, self.inner.clone(), None)
    }

    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        Self::readonly_positions(py, self.inner.clone(), None)
    }

    fn model_positions<'py>(
        &self,
        py: Python<'py>,
        model: PyModelIndex,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let (snapshot, local) = self
            .inner
            .model_snapshot(model.0)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("model index out of range"))?;
        Self::readonly_positions(py, snapshot, Some(local))
    }

    fn model_snapshot(&self, model: PyModelIndex) -> Option<(Self, PyModelIndex)> {
        self.inner
            .model_snapshot(model.0)
            .map(|(structure, local)| (Self::new(structure), PyModelIndex(local)))
    }

    fn ragged_models(&self) -> Option<Vec<Self>> {
        self.inner
            .ragged_models()
            .map(|models| models.iter().cloned().map(Self::new).collect())
    }

    #[getter]
    fn entity_count(&self) -> usize {
        self.inner.entity_count()
    }

    #[getter]
    fn generation(&self) -> PyCoordinateGeneration {
        PyCoordinateGeneration(self.inner.generation())
    }

    #[getter]
    fn annotations(&self) -> PyAtomAnnotations {
        PyAtomAnnotations(self.inner.annotations().clone())
    }

    #[getter]
    fn extensions(&self) -> PyExtensionStore {
        PyExtensionStore(self.inner.extensions().clone())
    }

    fn resolve(&self, symbol: PySymbolId) -> Option<String> {
        self.inner.resolve(symbol.0).map(str::to_owned)
    }

    fn without_extensions(&self) -> Self {
        Self::new(self.inner.without_extensions())
    }

    fn validate(&self) -> Vec<PyDiagnostic> {
        molframe::core::structure::validate(self.inner.data())
            .into_iter()
            .map(Into::into)
            .collect()
    }
}

#[cfg(test)]
#[path = "binding_tests.rs"]
mod tests;
