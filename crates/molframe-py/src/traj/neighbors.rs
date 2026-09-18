//! Reusable trajectory neighbour lists with a zero-copy `NumPy` query path.

use super::reader_types::PyTimestep;
use crate::crystallography::PyUnitCell;
use crate::geometry::borrowed_coordinates;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

type PyNeighborArrays<'py> = (
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<u32>>,
    Bound<'py, PyArray1<f32>>,
);

#[pyclass(name = "NeighborStatistics", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyNeighborStatistics {
    #[pyo3(get)]
    pub(crate) frames: usize,
    #[pyo3(get)]
    pub(crate) rebuilds: usize,
}

#[pyclass(name = "FrameNeighborList")]
pub(crate) struct PyFrameNeighborList {
    inner: molframe::traj::FrameNeighborList,
}

#[pymethods]
impl PyFrameNeighborList {
    #[new]
    fn new(left: Vec<u32>, right: Vec<u32>, cutoff: f32, skin: f32) -> Self {
        Self {
            inner: molframe::traj::FrameNeighborList::new(left, right, cutoff, skin),
        }
    }

    /// Compatibility entry point for one owned native timestep.
    ///
    /// Coordinates are borrowed from the timestep and the spatial kernel runs
    /// without the GIL. Prefer [`Self::pairs_positions`] for array workflows.
    fn pairs(&mut self, py: Python<'_>, timestep: &PyTimestep) -> PyResult<Vec<(u32, u32, f32)>> {
        let cell = timestep.cell.as_ref().map(|cell| cell.cell);
        py.detach(|| self.inner.pairs_positions(&timestep.positions, cell))
            .map(|pairs| {
                pairs
                    .into_iter()
                    .map(|pair| (pair.first, pair.second, pair.distance_squared))
                    .collect()
            })
            .map_err(neighbor_error)
    }

    /// Return three contiguous result arrays from C-contiguous `(atoms, 3)` input.
    ///
    /// The coordinate input is borrowed for this native call. The output tuple
    /// is `(first, second, distance_squared)` with dtypes `uint32`, `uint32`,
    /// and `float32`; it avoids a Python object per neighbour pair.
    #[pyo3(signature = (positions, *, cell=None))]
    fn pairs_positions<'py>(
        &mut self,
        py: Python<'py>,
        positions: PyReadonlyArray2<'_, f32>,
        cell: Option<&PyUnitCell>,
    ) -> PyResult<PyNeighborArrays<'py>> {
        let positions = borrowed_coordinates(&positions)?;
        let cell = cell.map(|cell| cell.cell);
        let pairs = py
            .detach(|| self.inner.pairs_positions(positions, cell))
            .map_err(neighbor_error)?;
        let mut first = Vec::with_capacity(pairs.len());
        let mut second = Vec::with_capacity(pairs.len());
        let mut distance_squared = Vec::with_capacity(pairs.len());
        for pair in pairs {
            first.push(pair.first);
            second.push(pair.second);
            distance_squared.push(pair.distance_squared);
        }
        Ok((
            first.into_pyarray(py),
            second.into_pyarray(py),
            distance_squared.into_pyarray(py),
        ))
    }

    #[getter]
    fn statistics(&self) -> PyNeighborStatistics {
        let value = self.inner.statistics();
        PyNeighborStatistics {
            frames: value.frames,
            rebuilds: value.rebuilds,
        }
    }
}

fn neighbor_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNeighborStatistics>()?;
    module.add_class::<PyFrameNeighborList>()
}
