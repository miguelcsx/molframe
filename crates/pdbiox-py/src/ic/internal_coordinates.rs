//! Python-owned projections for internal-coordinate measurement and rebuilding.

use crate::errors::read_error;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "Hedron", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHedron {
    #[pyo3(get)]
    first_length: f64,
    #[pyo3(get)]
    angle: f64,
    #[pyo3(get)]
    second_length: f64,
}

#[pymethods]
impl PyHedron {
    #[staticmethod]
    fn from_points(i: [f32; 3], j: [f32; 3], k: [f32; 3]) -> Option<Self> {
        pdbiox::Hedron::from_points(i, j, k).map(Self::from)
    }
}

#[pyclass(name = "Dihedron", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDihedron {
    #[pyo3(get)]
    first: PyHedron,
    #[pyo3(get)]
    angle: f64,
    #[pyo3(get)]
    length: f64,
    #[pyo3(get)]
    torsion: f64,
}

#[pymethods]
impl PyDihedron {
    #[staticmethod]
    fn from_points(i: [f32; 3], j: [f32; 3], k: [f32; 3], l: [f32; 3]) -> Option<Self> {
        pdbiox::Dihedron::from_points(i, j, k, l).map(Self::from)
    }
}

#[pyclass(name = "InternalAtom", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyInternalAtom {
    #[pyo3(get)]
    atom: u32,
    #[pyo3(get)]
    references: [u32; 3],
    #[pyo3(get)]
    coordinate: PyDihedron,
}

#[pyclass(name = "BatFrame", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBatFrame {
    inner: pdbiox::BatFrame,
}

#[pymethods]
impl PyBatFrame {
    #[getter]
    fn seed_positions(&self) -> Vec<[f32; 3]> {
        self.inner.seed_positions().to_vec()
    }

    #[getter]
    fn coordinates(&self) -> Vec<[f64; 3]> {
        self.inner.coordinates().to_vec()
    }
}

#[pyclass(name = "InternalCoordinates", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyInternalCoordinates {
    inner: pdbiox::InternalCoordinates,
}

#[pymethods]
impl PyInternalCoordinates {
    #[getter]
    fn seeds(&self) -> Vec<(u32, [f32; 3])> {
        self.inner
            .seeds()
            .iter()
            .map(|(atom, position)| (atom.get(), *position))
            .collect()
    }

    #[getter]
    fn atoms(&self) -> Vec<PyInternalAtom> {
        self.inner.atoms().iter().copied().map(Into::into).collect()
    }

    fn measure_bat(
        &self,
        py: Python<'_>,
        positions: Vec<Option<[f32; 3]>>,
    ) -> PyResult<PyBatFrame> {
        let inner = self.inner.clone();
        py.detach(move || inner.measure_bat(&positions))
            .map(|inner| PyBatFrame { inner })
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }

    fn rebuild_bat(&self, py: Python<'_>, frame: &PyBatFrame) -> PyResult<Vec<Option<[f32; 3]>>> {
        let inner = self.inner.clone();
        let frame = frame.inner.clone();
        py.detach(move || inner.rebuild_bat(&frame))
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }

    fn rebuild(&self, py: Python<'_>) -> PyResult<Vec<Option<[f32; 3]>>> {
        let inner = self.inner.clone();
        py.detach(move || inner.rebuild())
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }
}

#[pyfunction]
#[pyo3(signature = (structure, model=0))]
pub(crate) fn internal_coordinates(
    py: Python<'_>,
    structure: &PyStructure,
    model: u32,
) -> PyResult<PyInternalCoordinates> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::internal_coordinates(&structure, pdbiox::ModelIndex::new(model)))
        .map(|inner| PyInternalCoordinates { inner })
        .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
}

#[pyfunction]
pub(crate) fn place_atom(
    i: [f32; 3],
    j: [f32; 3],
    k: [f32; 3],
    length: f64,
    angle: f64,
    torsion: f64,
) -> Option<[f32; 3]> {
    pdbiox::place_atom(i, j, k, length, angle, torsion)
}

impl From<pdbiox::Hedron> for PyHedron {
    fn from(value: pdbiox::Hedron) -> Self {
        Self {
            first_length: value.first_length,
            angle: value.angle,
            second_length: value.second_length,
        }
    }
}

impl From<pdbiox::Dihedron> for PyDihedron {
    fn from(value: pdbiox::Dihedron) -> Self {
        Self {
            first: value.first.into(),
            angle: value.angle,
            length: value.length,
            torsion: value.torsion,
        }
    }
}

impl From<pdbiox::InternalAtom> for PyInternalAtom {
    fn from(value: pdbiox::InternalAtom) -> Self {
        Self {
            atom: value.atom.get(),
            references: value.references.map(pdbiox::AtomIndex::get),
            coordinate: value.coordinate.into(),
        }
    }
}

#[cfg(test)]
#[path = "internal_coordinates_tests.rs"]
mod tests;
