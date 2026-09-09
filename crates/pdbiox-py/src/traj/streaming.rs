//! Streaming numerical consumers with execution-owned result lifetimes.

use crate::core::execution::PyExecutionContext;
use numpy::{PyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pdbiox::core::ExecutionContext;
use pdbiox::traj::{FrameAlignment, TrajectoryReader, TrajectoryReaderOptions};
use pyo3::exceptions::{PyMemoryError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

fn open(
    path: &std::path::Path,
    context: &ExecutionContext,
    bytes: usize,
) -> PyResult<Box<dyn TrajectoryReader>> {
    pdbiox::traj::read_trajectory_in(
        path,
        &TrajectoryReaderOptions {
            memory_limit_bytes: bytes,
            ..TrajectoryReaderOptions::default()
        },
        context,
    )
    .map_err(|error| match error {
        pdbiox::traj::TrajectoryIoError::Reader(pdbiox::traj::TrajectoryError::MemoryLimit {
            ..
        }) => PyMemoryError::new_err(error.to_string()),
        _ => value_error(error),
    })
}

/// One-pass RMSF. The `NumPy` base retains its native allocation reservation.
#[pyfunction]
#[pyo3(signature = (path, context, *, frame_workspace_bytes=8_000_000))]
fn rmsf_stream<'py>(
    py: Python<'py>,
    path: PathBuf,
    context: &PyExecutionContext,
    frame_workspace_bytes: usize,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let context = context.native();
    let result = py.detach(move || {
        let mut reader = open(&path, &context, frame_workspace_bytes)?;
        pdbiox::traj::rmsf_stream(&mut *reader, &context, frame_workspace_bytes)
            .map_err(trajectory_error)
    })?;
    crate::core::retained_array::retained_scalars(py, result)
}

/// Emits RMSD rows incrementally through a Python callback.
#[pyfunction]
#[pyo3(signature = (path, reference, emit, context, *, fit=true, frame_workspace_bytes=8_000_000))]
fn rmsd_stream(
    py: Python<'_>,
    path: PathBuf,
    reference: PyReadonlyArray2<'_, f32>,
    emit: Py<PyAny>,
    context: &PyExecutionContext,
    fit: bool,
    frame_workspace_bytes: usize,
) -> PyResult<u64> {
    if reference.shape()[1] != 3 {
        return Err(PyValueError::new_err(
            "reference must have shape (atoms, 3)",
        ));
    }
    let (positions, _) = reference.as_slice()?.as_chunks::<3>();
    let context = context.native();
    py.detach(move || {
        let mut reader = open(&path, &context, frame_workspace_bytes)?;
        let mut callback_error = None;
        let result = pdbiox::traj::rmsd_stream(
            &mut *reader,
            positions,
            if fit {
                FrameAlignment::Rigid
            } else {
                FrameAlignment::None
            },
            &context,
            frame_workspace_bytes,
            |index, time, rmsd| {
                Python::attach(|py| emit.call1(py, (index, time, rmsd)))
                    .map(|_| ())
                    .map_err(|error| {
                        callback_error = Some(error);
                        pdbiox::traj::TrajectoryError::InvalidSource {
                            format: "Python RMSD sink",
                        }
                    })
            },
        );
        match callback_error {
            Some(error) => Err(error),
            None => result.map_err(trajectory_error),
        }
    })
}

/// Emits contact counts through the native bounded frame consumer.
#[pyfunction]
#[pyo3(signature = (path, cutoff, emit, context, *, options=None, frame_workspace_bytes=8_000_000))]
fn contact_counts_stream(
    py: Python<'_>,
    path: PathBuf,
    cutoff: f32,
    emit: Py<PyAny>,
    context: &PyExecutionContext,
    options: Option<crate::spatial::PySpatialSearchOptions>,
    frame_workspace_bytes: usize,
) -> PyResult<u64> {
    let context = context.native();
    let options = options.map_or_else(pdbiox::SpatialSearchOptions::default, |options| {
        options.inner()
    });
    py.detach(move || {
        let mut reader = open(&path, &context, frame_workspace_bytes)?;
        let mut callback_error = None;
        let result = pdbiox::traj::contact_counts_stream(
            &mut *reader,
            cutoff,
            options,
            &context,
            frame_workspace_bytes,
            |index, time, count| {
                Python::attach(|py| emit.call1(py, (index, time, count)))
                    .map(|_| ())
                    .map_err(|error| {
                        callback_error = Some(error);
                        pdbiox::traj::TrajectoryError::InvalidSource {
                            format: "Python contact-count sink",
                        }
                    })
            },
        );
        match callback_error {
            Some(error) => Err(error),
            None => result.map_err(trajectory_error),
        }
    })
}

/// Runs a callback over leased frames without materializing a trajectory.
#[pyfunction]
#[pyo3(signature = (path, consume, context, *, frame_workspace_bytes=8_000_000))]
fn run_analysis_stream(
    py: Python<'_>,
    path: PathBuf,
    consume: Py<PyAny>,
    context: &PyExecutionContext,
    frame_workspace_bytes: usize,
) -> PyResult<u64> {
    use pdbiox::core::{Backpressure, BatchDemand, BatchSource, ChunkId, DatasetId, LogicalRow};
    let context = context.native();
    py.detach(move || {
        let reader = open(&path, &context, frame_workspace_bytes)?;
        let mut source = pdbiox::traj::TrajectoryBatchSource::new(
            reader,
            DatasetId::new(0),
            ChunkId::new(0),
            LogicalRow::new(0),
        );
        let mut count = 0_u64;
        loop {
            if context.cancellation().is_cancelled() {
                return Err(PyValueError::new_err("trajectory execution cancelled"));
            }
            match source
                .next_batch(BatchDemand::new(1, frame_workspace_bytes), &context)
                .map_err(trajectory_error)?
            {
                Backpressure::Finished => return Ok(count),
                Backpressure::Pending => {
                    return Err(PyMemoryError::new_err(
                        "retained frames exhaust the execution budget",
                    ));
                }
                Backpressure::Ready(lease) => {
                    Python::attach(|py| {
                        let frame = Bound::new(
                            py,
                            super::stream_frame::PyStreamFrame {
                                lease: std::sync::Arc::new(lease),
                            },
                        )?;
                        consume.call1(py, (frame,)).map(|_| ())
                    })?;
                    count = count
                        .checked_add(1)
                        .ok_or_else(|| PyValueError::new_err("frame count overflow"))?;
                }
            }
        }
    })
}

fn trajectory_error(error: pdbiox::traj::TrajectoryError) -> PyErr {
    match error {
        pdbiox::traj::TrajectoryError::MemoryLimit { .. }
        | pdbiox::traj::TrajectoryError::Spatial(pdbiox::SpatialError::Memory(_)) => {
            PyMemoryError::new_err(error.to_string())
        }
        _ => value_error(error),
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<super::stream_frame::PyStreamFrame>()?;
    module.add_function(wrap_pyfunction!(run_analysis_stream, module)?)?;
    module.add_function(wrap_pyfunction!(rmsf_stream, module)?)?;
    module.add_function(wrap_pyfunction!(rmsd_stream, module)?)?;
    module.add_function(wrap_pyfunction!(contact_counts_stream, module)?)?;
    Ok(())
}
