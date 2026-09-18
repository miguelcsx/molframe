//! Python projections for the reusable trajectory primitives.
//!
//! The native reader traits are intentionally adapted to owned Python sources:
//! a reader owns its frame store, while every `read_next` operation returns one
//! typed frame.  All iteration, validation and transforms stay in Rust.

use super::PyTrajectory;
use super::reader_types::PyTimestep;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "Frame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFrame {
    #[pyo3(get)]
    pub(crate) positions: Vec<[f32; 3]>,
}

#[pymethods]
impl PyFrame {
    #[new]
    fn new(positions: Vec<[f32; 3]>) -> Self {
        Self { positions }
    }

    #[getter]
    fn atom_count(&self) -> usize {
        self.positions.len()
    }
}

impl From<PyFrame> for molframe::traj::Frame {
    fn from(value: PyFrame) -> Self {
        Self {
            positions: value.positions,
        }
    }
}

impl From<molframe::traj::Frame> for PyFrame {
    fn from(value: molframe::traj::Frame) -> Self {
        Self {
            positions: value.positions,
        }
    }
}

#[pyclass(name = "Units", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyUnits {
    #[pyo3(get)]
    pub(crate) length: &'static str,
    #[pyo3(get)]
    pub(crate) time: &'static str,
    #[pyo3(get)]
    pub(crate) force: &'static str,
}

impl From<molframe::traj::Units> for PyUnits {
    fn from(value: molframe::traj::Units) -> Self {
        Self {
            length: value.length,
            time: value.time,
            force: value.force,
        }
    }
}

#[pyclass(name = "RandomAccess", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyRandomAccess {
    Full,
    ViaIndex,
    None,
}

#[pyclass(name = "MinimalTopology", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMinimalTopology {
    atom_count: usize,
}

#[pymethods]
impl PyMinimalTopology {
    #[new]
    fn new(atom_count: usize) -> Self {
        Self { atom_count }
    }

    #[getter]
    const fn atom_count(&self) -> usize {
        self.atom_count
    }

    fn validate_reader(&self, reader: &PyMemoryReader) -> PyResult<()> {
        if self.atom_count != reader.atom_count {
            return Err(PyValueError::new_err(format!(
                "trajectory atom-count mismatch: expected {}, found {}",
                self.atom_count, reader.atom_count
            )));
        }
        Ok(())
    }
}

#[pyclass(name = "MemoryReader", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMemoryReader {
    pub(crate) frames: Vec<PyTimestep>,
    pub(crate) cursor: usize,
    pub(crate) atom_count: usize,
}

impl PyMemoryReader {
    fn build(frames: Vec<PyTimestep>) -> PyResult<Self> {
        let atom_count = frames.first().map_or(0, |frame| frame.positions.len());
        if frames
            .iter()
            .any(|frame| frame.positions.len() != atom_count)
        {
            return Err(PyValueError::new_err(
                "trajectory frames must have one atom count",
            ));
        }
        Ok(Self {
            frames,
            cursor: 0,
            atom_count,
        })
    }

    pub(crate) fn read_next_frame(&mut self) -> Option<PyTimestep> {
        let frame = self.frames.get(self.cursor).cloned()?;
        self.cursor += 1;
        Some(frame)
    }
}

#[pymethods]
impl PyMemoryReader {
    #[new]
    fn new(py: Python<'_>, trajectory: PyRef<'_, PyTrajectory>) -> PyResult<Self> {
        Self::build(trajectory.timesteps(py)?)
    }

    #[staticmethod]
    fn from_frames(frames: Vec<PyTimestep>) -> PyResult<Self> {
        Self::build(frames)
    }

    #[getter]
    const fn n_atoms(&self) -> usize {
        self.atom_count
    }

    #[getter]
    fn n_frames(&self) -> usize {
        self.frames.len()
    }

    #[getter]
    fn units(&self) -> PyUnits {
        molframe::traj::Units::CANONICAL.into()
    }

    #[getter]
    fn random_access(&self) -> PyRandomAccess {
        PyRandomAccess::Full
    }

    fn read_next(&mut self) -> Option<PyTimestep> {
        self.read_next_frame()
    }

    fn read_all(&mut self) -> Vec<PyTimestep> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_next_frame() {
            frames.push(frame);
        }
        frames
    }

    fn seek(&mut self, frame: usize) {
        self.cursor = frame.min(self.frames.len());
    }

    fn rewind(&mut self) {
        self.cursor = 0;
    }
}

#[pyclass(name = "StreamingReader", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStreamingReader {
    frames: Vec<PyTimestep>,
    cursor: usize,
    atom_count: usize,
}

#[pymethods]
impl PyStreamingReader {
    #[new]
    #[pyo3(signature = (frames, *, atom_count=None))]
    fn new(frames: Vec<PyFrame>, atom_count: Option<usize>) -> PyResult<Self> {
        let frames = frames
            .into_iter()
            .enumerate()
            .map(|(frame, value)| PyTimestep {
                frame,
                time: None,
                dt: None,
                positions: value.positions,
                velocities: None,
                forces: None,
                cell: None,
                data: Vec::new(),
            })
            .collect::<Vec<_>>();
        let actual = frames.first().map_or(0, |frame| frame.positions.len());
        let atom_count = match atom_count {
            Some(value) => value,
            None => actual,
        };
        if frames
            .iter()
            .any(|frame| frame.positions.len() != atom_count)
        {
            return Err(PyValueError::new_err(
                "stream frames must match the declared atom count",
            ));
        }
        Ok(Self {
            frames,
            cursor: 0,
            atom_count,
        })
    }

    #[getter]
    const fn n_atoms(&self) -> usize {
        self.atom_count
    }

    #[getter]
    const fn n_frames(&self) -> Option<usize> {
        None
    }

    #[getter]
    fn units(&self) -> PyUnits {
        molframe::traj::Units::CANONICAL.into()
    }

    #[getter]
    const fn random_access(&self) -> PyRandomAccess {
        PyRandomAccess::None
    }

    fn read_next(&mut self) -> Option<PyTimestep> {
        let frame = self.frames.get(self.cursor).cloned()?;
        self.cursor += 1;
        Some(frame)
    }

    fn seek(&mut self, frame: usize) -> PyResult<()> {
        let _ = frame;
        Err(PyValueError::new_err(
            "streaming reader does not support random access",
        ))
    }

    fn rewind(&mut self) -> PyResult<()> {
        self.seek(0)
    }
}

#[pyclass(name = "ChainedReader", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainedReader {
    readers: Vec<PyMemoryReader>,
    current: usize,
    frame: usize,
    atom_count: usize,
}

#[pymethods]
impl PyChainedReader {
    #[new]
    fn new(readers: Vec<PyMemoryReader>) -> PyResult<Self> {
        let atom_count = readers.first().map_or(0, |reader| reader.atom_count);
        if readers.iter().any(|reader| reader.atom_count != atom_count) {
            return Err(PyValueError::new_err(
                "chained readers must have one atom count",
            ));
        }
        Ok(Self {
            readers,
            current: 0,
            frame: 0,
            atom_count,
        })
    }

    #[getter]
    const fn n_atoms(&self) -> usize {
        self.atom_count
    }

    #[getter]
    fn n_frames(&self) -> Option<usize> {
        self.readers.iter().try_fold(0usize, |total, reader| {
            total.checked_add(reader.frames.len())
        })
    }

    #[getter]
    fn random_access(&self) -> PyRandomAccess {
        PyRandomAccess::None
    }

    fn read_next(&mut self) -> Option<PyTimestep> {
        loop {
            let reader = self.readers.get_mut(self.current)?;
            if let Some(mut frame) = reader.read_next_frame() {
                frame.frame = self.frame;
                self.frame += 1;
                return Some(frame);
            }
            self.current += 1;
        }
    }

    fn seek(&mut self, frame: usize) -> PyResult<()> {
        let _ = frame;
        Err(PyValueError::new_err(
            "chained reader does not support random access",
        ))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFrame>()?;
    module.add_class::<PyUnits>()?;
    module.add_class::<PyRandomAccess>()?;
    module.add_class::<PyMinimalTopology>()?;
    module.add_class::<PyMemoryReader>()?;
    module.add_class::<PyStreamingReader>()?;
    module.add_class::<PyChainedReader>()?;
    Ok(())
}
