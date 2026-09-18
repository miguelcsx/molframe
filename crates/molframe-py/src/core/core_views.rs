//! Zero-copy selection views with explicit gather semantics for sparse positions.

use crate::core_values::{PyAabb, PyCoordinateGeneration};
use crate::index::PyModelIndex;
use crate::query::PySelection;
use crate::structure::PyStructure;
use numpy::IntoPyArray;
use numpy::ndarray::Array2;
use pyo3::exceptions::PyOverflowError;
use pyo3::prelude::*;

#[pyclass(name = "StructureView", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructureView(pub(crate) molframe::StructureView);

#[pymethods]
impl PyStructureView {
    #[staticmethod]
    fn new(structure: &PyStructure, selection: &PySelection) -> Self {
        Self(molframe::StructureView::new(
            structure.structure(),
            selection.inner.clone(),
        ))
    }

    #[getter]
    fn selection(&self) -> PySelection {
        PySelection {
            inner: self.0.selection().clone(),
        }
    }

    fn __len__(&self) -> PyResult<usize> {
        usize::try_from(self.0.len())
            .map_err(|_| PyOverflowError::new_err("view length exceeds Python limits"))
    }

    fn len(&self) -> u64 {
        self.0.len()
    }

    fn data(&self) -> crate::core_data::PyStructureData {
        crate::core_data::PyStructureData(self.0.data().clone())
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[getter]
    fn generation(&self) -> PyCoordinateGeneration {
        PyCoordinateGeneration(self.0.generation())
    }

    fn is_stale_for(&self, structure: &PyStructure) -> bool {
        self.0.is_stale_for(structure.structure())
    }

    fn narrow(&self, selection: &PySelection) -> Self {
        Self(self.0.narrow(&selection.inner))
    }

    fn union(&self, other: &Self) -> Self {
        Self(self.0.union(&other.0))
    }

    fn intersect(&self, other: &Self) -> Self {
        Self(self.0.intersect(&other.0))
    }

    fn difference(&self, other: &Self) -> Self {
        Self(self.0.difference(&other.0))
    }

    fn positions<'py>(
        &self,
        py: Python<'py>,
        model: PyModelIndex,
    ) -> Bound<'py, numpy::PyArray2<f32>> {
        let view = self.0.clone();
        let values = py.detach(|| view.positions(model.0).flatten().collect::<Vec<_>>());
        let shape = (values.len() / 3, 3);
        match Array2::from_shape_vec(shape, values) {
            Ok(array) => array.into_pyarray(py),
            Err(_) => Array2::zeros((0, 3)).into_pyarray(py),
        }
    }

    fn bounds(&self, model: PyModelIndex) -> PyAabb {
        PyAabb(self.0.bounds(model.0))
    }
}

#[pymethods]
impl PyStructure {
    fn view(&self) -> PyStructureView {
        PyStructureView(self.structure().view())
    }

    fn view_of(&self, selection: &PySelection) -> PyStructureView {
        PyStructureView(self.structure().view_of(selection.inner.clone()))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStructureView>()?;
    Ok(())
}
