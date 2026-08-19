//! Immutable rigid transforms over arrays and structure snapshots.

use super::arrays::borrowed_coordinates;
use crate::errors::read_error;
use crate::query::PySelection;
use crate::structure::PyStructure;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;

#[pyclass(name = "Rigid", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyRigid(pub(crate) pdbiox::Rigid);

#[pymethods]
impl PyRigid {
    #[classattr]
    const IDENTITY: Self = Self(pdbiox::Rigid::IDENTITY);

    #[new]
    fn new(rotation: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self(pdbiox::Rigid::new(rotation, translation))
    }

    #[staticmethod]
    fn identity() -> Self {
        Self(pdbiox::Rigid::IDENTITY)
    }

    #[staticmethod]
    fn translation(offset: [f64; 3]) -> Self {
        Self(pdbiox::Rigid::translation(offset))
    }

    #[getter]
    fn rotation(&self) -> [[f64; 3]; 3] {
        self.0.rotation
    }

    #[getter]
    fn offset(&self) -> [f64; 3] {
        self.0.translation
    }

    fn inverse(&self) -> Self {
        Self(self.0.inverse())
    }

    fn then(&self, next: &Self) -> Self {
        Self(self.0.then(&next.0))
    }

    fn apply<'py>(
        &self,
        py: Python<'py>,
        positions: PyReadonlyArray2<'_, f32>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let positions = borrowed_coordinates(&positions)?;
        let transform = self.0;
        let values = py.detach(|| {
            positions
                .iter()
                .flat_map(|&position| transform.apply(position))
                .collect::<Vec<_>>()
        });
        Array2::from_shape_vec((positions.len(), 3), values)
            .map(|array| array.into_pyarray(py))
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
    }

    fn apply_one(&self, position: [f32; 3]) -> [f32; 3] {
        self.0.apply(position)
    }

    fn apply_all<'py>(
        &self,
        py: Python<'py>,
        positions: PyReadonlyArray2<'_, f32>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let positions = borrowed_coordinates(&positions)?;
        let transform = self.0;
        let values = py.detach(|| {
            positions
                .iter()
                .flat_map(|&position| transform.apply(position))
                .collect::<Vec<_>>()
        });
        Array2::from_shape_vec((positions.len(), 3), values)
            .map(|array| array.into_pyarray(py))
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
    }

    fn determinant(&self) -> f64 {
        self.0.determinant()
    }
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (transform, *, selection=None))]
    fn transformed(
        &self,
        py: Python<'_>,
        transform: &PyRigid,
        selection: Option<&PySelection>,
    ) -> PyResult<Self> {
        let structure = self.structure().clone();
        let selection = selection.map_or_else(
            || pdbiox::AtomSelection::All(structure.atom_count()),
            |selection| selection.inner.clone(),
        );
        let transform = transform.0;
        py.detach(|| pdbiox::transform(&structure, &selection, &transform))
            .map(Self::new)
            .map_err(|findings| read_error(py, &findings))
    }
}

#[pyclass(name = "Superposition", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySuperposition {
    pub(super) transform: PyRigid,
    pub(super) rmsd: f64,
}

#[pymethods]
impl PySuperposition {
    #[getter]
    const fn transform(&self) -> PyRigid {
        self.transform
    }
    #[getter]
    const fn rmsd(&self) -> f64 {
        self.rmsd
    }
}
