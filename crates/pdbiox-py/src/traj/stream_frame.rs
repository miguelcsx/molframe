//! Python ownership of a budgeted trajectory frame and its exported arrays.

use numpy::PyArray2;
use pdbiox::core::BatchLease;
use pdbiox::traj::{Timestep, TrajectoryBatch};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyclass(name = "StreamFrame", frozen, skip_from_py_object)]
pub(super) struct PyStreamFrame {
    pub(super) lease: Arc<BatchLease<TrajectoryBatch>>,
}

impl PyStreamFrame {
    fn frame(&self) -> PyResult<&Timestep> {
        self.lease
            .batch()
            .timestep()
            .ok_or_else(|| PyRuntimeError::new_err("frame is absent"))
    }
}

#[pymethods]
impl PyStreamFrame {
    #[getter(frame)]
    fn index(&self) -> PyResult<usize> {
        Ok(self.frame()?.frame)
    }
    #[getter]
    fn time(&self) -> PyResult<Option<f64>> {
        Ok(self.frame()?.time)
    }
    #[getter]
    fn dt(&self) -> PyResult<Option<f64>> {
        Ok(self.frame()?.dt)
    }
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        super::arrays::frame_vectors(py, Arc::clone(&self.lease), 0)?
            .ok_or_else(|| PyRuntimeError::new_err("positions are absent"))
    }
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f32>>>> {
        super::arrays::frame_vectors(py, Arc::clone(&self.lease), 1)
    }
    #[getter]
    fn forces<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f32>>>> {
        super::arrays::frame_vectors(py, Arc::clone(&self.lease), 2)
    }
    #[getter]
    fn cell(&self) -> PyResult<Option<crate::crystallography::PyUnitCell>> {
        self.frame()?
            .cell
            .map(crate::crystallography::PyUnitCell::from_native)
            .transpose()
    }
}
