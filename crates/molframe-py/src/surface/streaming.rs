//! Thin Python consumers for charged native surface streams.

use crate::core::execution::PyExecutionContext;
use crate::crystallography::PyUnitCell;
use crate::geometry::borrowed_coordinates;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyMemoryError;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyfunction]
#[pyo3(signature = (positions, radii, probe, samples, context, *, cell=None))]
fn collect_shrake_rupley<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    samples: u16,
    context: &PyExecutionContext,
    cell: Option<&PyUnitCell>,
) -> PyResult<Bound<'py, numpy::PyArray1<f64>>> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    let periodic = cell
        .map(|cell| molframe::PeriodicBox::from_cell(cell.cell))
        .transpose()
        .map_err(|error| surface_error(error.into()))?;
    let context = context.native();
    let result = py
        .detach(move || {
            molframe::surface::collect_shrake_rupley(
                positions,
                radii,
                probe,
                samples,
                periodic.as_ref(),
                &context,
            )
        })
        .map_err(surface_error)?;
    crate::core::retained_array::retained_scalars(py, result)
}

#[pyfunction]
#[pyo3(signature = (positions, radii, probe, samples, emit, context, *, cell=None))]
fn visit_shrake_rupley(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    samples: u16,
    emit: Py<PyAny>,
    context: &PyExecutionContext,
    cell: Option<&PyUnitCell>,
) -> PyResult<()> {
    let positions = borrowed_coordinates(&positions)?;
    let radii = radii.as_slice()?;
    let periodic = cell
        .map(|cell| molframe::PeriodicBox::from_cell(cell.cell))
        .transpose()
        .map_err(|error| surface_error(error.into()))?;
    let context = context.native();
    py.detach(move || {
        let mut callback_error = None;
        let result = molframe::surface::visit_shrake_rupley(
            positions,
            radii,
            probe,
            samples,
            periodic.as_ref(),
            &context,
            |atom, area| {
                Python::attach(|py| emit.call1(py, (atom, area)))
                    .map(|_| ())
                    .map_err(|error| {
                        callback_error = Some(error);
                        molframe::surface::SasaError::from(molframe::SpatialError::Cancelled)
                    })
            },
        );
        match callback_error {
            Some(error) => Err(error),
            None => result.map_err(surface_error),
        }
    })
}

#[pyfunction]
#[pyo3(signature = (path, radii, probe, samples, emit, context, *, frame_workspace_bytes=8_000_000))]
fn sasa_stream(
    py: Python<'_>,
    path: PathBuf,
    radii: PyReadonlyArray1<'_, f32>,
    probe: f32,
    samples: u16,
    emit: Py<PyAny>,
    context: &PyExecutionContext,
    frame_workspace_bytes: usize,
) -> PyResult<u64> {
    let radii = radii.as_slice()?;
    let context = context.native();
    py.detach(move || {
        let mut reader = molframe::traj::read_trajectory_in(
            &path,
            &molframe::traj::TrajectoryReaderOptions {
                memory_limit_bytes: frame_workspace_bytes,
                ..molframe::traj::TrajectoryReaderOptions::default()
            },
            &context,
        )
        .map_err(|error| match error {
            molframe::traj::TrajectoryIoError::Reader(
                molframe::traj::TrajectoryError::MemoryLimit { .. },
            ) => PyMemoryError::new_err(error.to_string()),
            _ => crate::core::errors::SasaStreamError::new_err(error.to_string()),
        })?;
        let mut callback_error = None;
        let result = molframe::analysis::sasa_stream(
            &mut *reader,
            radii,
            probe,
            samples,
            &context,
            frame_workspace_bytes,
            |frame, time, area| {
                Python::attach(|py| emit.call1(py, (frame, time, area)))
                    .map(|_| ())
                    .map_err(|error| {
                        callback_error = Some(error);
                        molframe::traj::TrajectoryError::Cancelled.into()
                    })
            },
        );
        match callback_error {
            Some(error) => Err(error),
            None => result.map_err(stream_error),
        }
    })
}

fn surface_error(error: molframe::surface::SasaError) -> PyErr {
    match error {
        molframe::surface::SasaError::Spatial(molframe::SpatialError::Memory(_)) => {
            PyMemoryError::new_err(error.to_string())
        }
        _ => crate::core::errors::SasaError::new_err(error.to_string()),
    }
}

fn stream_error(error: molframe::analysis::SasaStreamError) -> PyErr {
    match error {
        molframe::analysis::SasaStreamError::Surface(molframe::surface::SasaError::Spatial(
            molframe::SpatialError::Memory(_),
        ))
        | molframe::analysis::SasaStreamError::Trajectory(
            molframe::traj::TrajectoryError::MemoryLimit { .. },
        ) => PyMemoryError::new_err(error.to_string()),
        _ => crate::core::errors::SasaStreamError::new_err(error.to_string()),
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(visit_shrake_rupley, module)?)?;
    module.add_function(wrap_pyfunction!(collect_shrake_rupley, module)?)?;
    module.add_function(wrap_pyfunction!(sasa_stream, module)?)
}
