//! Streaming trajectory analyses: constant memory in the number of frames.
//!
//! Each function pulls frames through one reusable buffer, with the GIL released
//! and the run interruptible, and keeps only the per-frame results it returns.

use crate::trajectory_stream::Source;
use molframe::spatial::SpatialSearchOptions;
use molframe::trajectory::{
    DenseChunk, FrameAlignment, Timestep, TrajectoryError, contact_counts_stream, read_dense_chunk,
    rmsd_stream, rmsf_stream, run_analysis_stream,
};
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PySlice;
use std::path::PathBuf;

/// One value per streamed frame, with the frame it belongs to.
#[derive(Debug)]
#[pyclass(
    name = "FrameSeries",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
struct PyFrameSeries {
    frames: Py<PyArray1<u64>>,
    times: Py<PyArray1<f64>>,
    values: Py<PyAny>,
}

#[pymethods]
impl PyFrameSeries {
    /// Source frame ordinal of every value.
    #[getter]
    fn frames(&self, py: Python<'_>) -> Py<PyArray1<u64>> {
        self.frames.clone_ref(py)
    }

    /// Simulation time of every value in picoseconds; `nan` where unrecorded.
    #[getter]
    fn times(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.times.clone_ref(py)
    }

    /// The analysis value of every frame.
    #[getter]
    fn values(&self, py: Python<'_>) -> Py<PyAny> {
        self.values.clone_ref(py)
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.frames.bind(py).len()
    }
}

/// Per-frame results gathered while streaming.
#[derive(Default)]
struct Collected<T> {
    frames: Vec<u64>,
    times: Vec<f64>,
    values: Vec<T>,
}

impl<T> Collected<T> {
    fn push(&mut self, frame: usize, time: Option<f64>, value: T) -> Result<(), TrajectoryError> {
        self.frames
            .push(u64::try_from(frame).map_err(|_| TrajectoryError::IdentityOverflow)?);
        self.times.push(time.unwrap_or(f64::NAN));
        self.values.push(value);
        Ok(())
    }

    fn into_series(
        self,
        py: Python<'_>,
        values: impl FnOnce(Python<'_>, Vec<T>) -> Py<PyAny>,
    ) -> PyFrameSeries {
        let frames = PyArray1::from_vec(py, self.frames);
        let times = PyArray1::from_vec(py, self.times);
        frames.readwrite().make_nonwriteable();
        times.readwrite().make_nonwriteable();
        PyFrameSeries {
            frames: frames.unbind(),
            times: times.unbind(),
            values: values(py, self.values),
        }
    }
}

fn readonly_values<T: numpy::Element>(py: Python<'_>, values: Vec<T>) -> Py<PyAny> {
    let array = PyArray1::from_vec(py, values);
    array.readwrite().make_nonwriteable();
    array.into_any().unbind()
}

/// The reference positions: an `(atoms, 3)` array, or a source frame index read
/// through the same atom subset.
fn reference_positions(
    source: &Source,
    reference: &Bound<'_, PyAny>,
    atoms: Option<&Vec<usize>>,
) -> PyResult<Vec<[f32; 3]>> {
    let py = reference.py();
    if let Ok(index) = reference.extract::<usize>() {
        let out_of_range = || crate::error::index("reference frame out of range");
        let start = isize::try_from(index).map_err(|_| out_of_range())?;
        let slice = PySlice::new(py, start, start.checked_add(1).ok_or_else(out_of_range)?, 1);
        let (mut reader, _) = source.open(Some(&slice), atoms.cloned())?;
        let chunk: DenseChunk = read_dense_chunk(&mut reader, &mut Timestep::default(), 1)
            .map_err(crate::error::kernel)?;
        if chunk.frames() == 0 {
            return Err(crate::error::index("reference frame out of range"));
        }
        return Ok(chunk.coordinates.as_chunks::<3>().0.to_vec());
    }
    let array = py
        .import("numpy")?
        .call_method1("ascontiguousarray", (reference, "float32"))?
        .cast_into::<PyArray2<f32>>()
        .map_err(|_| crate::error::value("reference must have shape (atoms, 3)"))?;
    let readonly = array.readonly();
    if readonly.shape()[1] != 3 {
        return Err(crate::error::value("reference must have shape (atoms, 3)"));
    }
    let flat = readonly
        .as_slice()
        .map_err(|_| crate::error::value("reference must be C-contiguous"))?;
    Ok(flat.as_chunks::<3>().0.to_vec())
}

/// RMSD of every streamed frame against one reference, in constant memory.
///
/// `reference` is an `(atoms, 3)` array (any dtype or strides) or a source frame
/// index. With `align=True` each frame is rigidly fitted onto it first. `frames`,
/// `atoms` and `memory_limit` mean what they do for `read`.
#[pyfunction]
#[pyo3(signature = (path, reference, *, align=true, format=None, frames=None, atoms=None, memory_limit=None, context=None))]
#[allow(clippy::too_many_arguments)]
#[pyo3(name = "rmsd_stream")]
fn rmsd_stream_py(
    py: Python<'_>,
    path: PathBuf,
    reference: &Bound<'_, PyAny>,
    align: bool,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<Vec<usize>>,
    memory_limit: Option<usize>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyFrameSeries> {
    let source = Source::new(path, format, memory_limit)?;
    let reference = reference_positions(&source, reference, atoms.as_ref())?;
    let (reader, bytes) = source.open(frames, atoms)?;
    let alignment = if align {
        FrameAlignment::Rigid
    } else {
        FrameAlignment::None
    };
    let collected = crate::execution::run(py, context, move |context| {
        let mut reader = reader;
        let mut collected = Collected::default();
        rmsd_stream(
            &mut reader,
            &reference,
            alignment,
            context,
            bytes,
            |frame, time, value| collected.push(frame, time, value),
        )?;
        Ok::<_, TrajectoryError>(collected)
    })?
    .map_err(crate::error::kernel)?;
    Ok(collected.into_series(py, readonly_values))
}

/// Per-atom RMSF over the streamed frames, in constant memory.
///
/// No alignment is applied: fit the trajectory first if it drifts or rotates.
#[pyfunction]
#[pyo3(signature = (path, *, format=None, frames=None, atoms=None, memory_limit=None, context=None))]
#[pyo3(name = "rmsf_stream")]
fn rmsf_stream_py<'py>(
    py: Python<'py>,
    path: PathBuf,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<Vec<usize>>,
    memory_limit: Option<usize>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let source = Source::new(path, format, memory_limit)?;
    let (reader, bytes) = source.open(frames, atoms)?;
    let values = crate::execution::run(py, context, move |context| {
        let mut reader = reader;
        rmsf_stream(&mut reader, context, bytes).map(|values| values.to_vec())
    })?
    .map_err(crate::error::kernel)?;
    let array = PyArray1::from_vec(py, values);
    array.readwrite().make_nonwriteable();
    Ok(array)
}

/// Unique atom pairs within `cutoff` for every streamed frame, in constant memory.
///
/// A frame's periodic cell is honoured when it has one. `backend` picks the
/// spatial search (`"auto"`, `"grid"`, ...).
#[pyfunction]
#[pyo3(signature = (path, cutoff, *, backend="auto", format=None, frames=None, atoms=None, memory_limit=None, context=None))]
#[allow(clippy::too_many_arguments)]
#[pyo3(name = "contact_counts_stream")]
fn contact_counts_stream_py(
    py: Python<'_>,
    path: PathBuf,
    cutoff: f32,
    backend: &str,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<Vec<usize>>,
    memory_limit: Option<usize>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyFrameSeries> {
    let options = SpatialSearchOptions {
        backend: crate::backend::parse(backend)?,
        ..SpatialSearchOptions::default()
    };
    let source = Source::new(path, format, memory_limit)?;
    let (reader, bytes) = source.open(frames, atoms)?;
    let collected = crate::execution::run(py, context, move |context| {
        let mut reader = reader;
        let mut collected = Collected::default();
        contact_counts_stream(
            &mut reader,
            cutoff,
            options,
            context,
            bytes,
            |frame, time, count| collected.push(frame, time, count),
        )?;
        Ok::<_, TrajectoryError>(collected)
    })?
    .map_err(crate::error::kernel)?;
    Ok(collected.into_series(py, readonly_values))
}

/// Why a callback-driven stream stopped.
enum Stop {
    Reader(TrajectoryError),
    Callback(PyErr),
}

impl From<TrajectoryError> for Stop {
    fn from(error: TrajectoryError) -> Self {
        Self::Reader(error)
    }
}

/// Calls `callback(frame, time, positions)` for every streamed frame; returns the count.
///
/// `positions` is a read-only `(atoms, 3)` float32 copy that outlives the call. The
/// decode runs with the GIL released and takes it only around each callback, so
/// memory stays at one frame however long the trajectory is. An exception in the
/// callback stops the stream and propagates.
#[pyfunction]
#[pyo3(signature = (path, callback, *, format=None, frames=None, atoms=None, memory_limit=None, context=None))]
#[pyo3(name = "run_analysis_stream")]
fn run_analysis_stream_py(
    py: Python<'_>,
    path: PathBuf,
    callback: Py<PyAny>,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<Vec<usize>>,
    memory_limit: Option<usize>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<u64> {
    let source = Source::new(path, format, memory_limit)?;
    let (reader, bytes) = source.open(frames, atoms)?;
    let outcome = crate::execution::run(py, context, move |context| {
        let mut reader = reader;
        run_analysis_stream(&mut reader, context, bytes, |frame, _| {
            Python::attach(|py| -> PyResult<()> {
                let positions = PyArray1::from_slice(py, frame.positions.as_flattened())
                    .reshape([frame.positions.len(), 3])?;
                positions.readwrite().make_nonwriteable();
                callback.call1(py, (frame.frame, frame.time, positions))?;
                Ok(())
            })
            .map_err(Stop::Callback)
        })
    })?;
    match outcome {
        Ok(count) => Ok(count),
        Err(Stop::Callback(error)) => Err(error),
        Err(Stop::Reader(error)) => Err(crate::error::kernel(error)),
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFrameSeries>()?;
    module.add_function(wrap_pyfunction!(rmsd_stream_py, module)?)?;
    module.add_function(wrap_pyfunction!(rmsf_stream_py, module)?)?;
    module.add_function(wrap_pyfunction!(contact_counts_stream_py, module)?)?;
    module.add_function(wrap_pyfunction!(run_analysis_stream_py, module)?)
}
