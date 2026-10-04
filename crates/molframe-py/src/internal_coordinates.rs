//! Mechanical bindings for the native internal-coordinate forest.

use crate::bindings::PyStructure;
use molframe::ic::{BatFrame, InternalCoordinates};
use pyo3::prelude::*;

#[derive(Clone, Debug)]
#[pyclass(
    name = "BatFrame",
    frozen,
    skip_from_py_object,
    module = "molframe.geometry"
)]
pub(crate) struct PyBatFrame {
    inner: BatFrame,
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

#[derive(Clone, Debug)]
#[pyclass(
    name = "InternalCoordinates",
    frozen,
    skip_from_py_object,
    module = "molframe.geometry"
)]
pub(crate) struct PyInternalCoordinates {
    inner: InternalCoordinates,
}

#[pymethods]
impl PyInternalCoordinates {
    #[getter]
    fn seeds(&self) -> Vec<(usize, [f32; 3])> {
        self.inner
            .seeds()
            .iter()
            .map(|(atom, position)| (atom.as_usize(), *position))
            .collect()
    }

    #[getter]
    fn atoms(&self) -> Vec<(usize, [usize; 3], [f64; 6])> {
        self.inner
            .atoms()
            .iter()
            .map(|atom| {
                let coordinate = atom.coordinate;
                (
                    atom.atom.as_usize(),
                    atom.references.map(molframe::AtomIndex::as_usize),
                    [
                        coordinate.first.first_length,
                        coordinate.first.angle,
                        coordinate.first.second_length,
                        coordinate.angle,
                        coordinate.length,
                        coordinate.torsion,
                    ],
                )
            })
            .collect()
    }

    fn rebuild(&self, py: Python<'_>) -> PyResult<Vec<Option<[f32; 3]>>> {
        py.detach(|| self.inner.rebuild())
            .map_err(crate::error::kernel)
    }

    fn measure_bat(
        &self,
        py: Python<'_>,
        positions: Vec<Option<[f32; 3]>>,
    ) -> PyResult<PyBatFrame> {
        let positions = positions.into_boxed_slice();
        py.detach(|| self.inner.measure_bat(&positions))
            .map(|inner| PyBatFrame { inner })
            .map_err(crate::error::kernel)
    }

    fn rebuild_bat(&self, py: Python<'_>, frame: &PyBatFrame) -> PyResult<Vec<Option<[f32; 3]>>> {
        py.detach(|| self.inner.rebuild_bat(&frame.inner))
            .map_err(crate::error::kernel)
    }
}

#[pyfunction]
#[pyo3(signature = (structure, *, model=0))]
fn internal_coordinates(
    py: Python<'_>,
    structure: &PyStructure,
    model: u32,
) -> PyResult<PyInternalCoordinates> {
    py.detach(|| {
        molframe::ic::internal_coordinates(
            structure.inner.engine(),
            molframe::ModelIndex::new(model),
        )
    })
    .map(|inner| PyInternalCoordinates { inner })
    .map_err(crate::error::kernel)
}

#[pyfunction]
fn place_atom(
    py: Python<'_>,
    first: [f32; 3],
    second: [f32; 3],
    third: [f32; 3],
    length: f64,
    angle: f64,
    torsion: f64,
) -> Option<[f32; 3]> {
    py.detach(|| molframe::ic::place_atom(first, second, third, length, angle, torsion))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyBatFrame>()?;
    module.add_class::<PyInternalCoordinates>()?;
    module.add_function(wrap_pyfunction!(internal_coordinates, module)?)?;
    module.add_function(wrap_pyfunction!(place_atom, module)?)?;
    Ok(())
}
