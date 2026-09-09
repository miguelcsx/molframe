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
        .map(|cell| pdbiox::PeriodicBox::from_cell(cell.cell))
        .transpose()
        .map_err(|error| surface_error(error.into()))?;
    let context = context.native();
    let result = py
        .detach(move || {
            pdbiox::surface::collect_shrake_rupley(
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
        .map(|cell| pdbiox::PeriodicBox::from_cell(cell.cell))
        .transpose()
        .map_err(|error| surface_error(error.into()))?;
    let context = context.native();
    py.detach(move || {
        let mut callback_error = None;
        let result = pdbiox::surface::visit_shrake_rupley(
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
                        pdbiox::surface::SasaError::from(pdbiox::SpatialError::Cancelled)
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
        let mut reader = pdbiox::traj::read_trajectory_in(
            &path,
            &pdbiox::traj::TrajectoryReaderOptions {
                memory_limit_bytes: frame_workspace_bytes,
                ..pdbiox::traj::TrajectoryReaderOptions::default()
            },
            &context,
        )
        .map_err(|error| match error {
            pdbiox::traj::TrajectoryIoError::Reader(
                pdbiox::traj::TrajectoryError::MemoryLimit { .. },
            ) => PyMemoryError::new_err(error.to_string()),
            _ => crate::core::errors::SasaStreamError::new_err(error.to_string()),
        })?;
        let mut callback_error = None;
        let result = pdbiox::analysis::sasa_stream(
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
                        pdbiox::traj::TrajectoryError::Cancelled.into()
                    })
            },
        );
        match callback_error {
            Some(error) => Err(error),
            None => result.map_err(stream_error),
        }
    })
}

fn surface_error(error: pdbiox::surface::SasaError) -> PyErr {
    match error {
        pdbiox::surface::SasaError::Spatial(pdbiox::SpatialError::Memory(_)) => {
            PyMemoryError::new_err(error.to_string())
        }
        _ => crate::core::errors::SasaError::new_err(error.to_string()),
    }
}

fn stream_error(error: pdbiox::analysis::SasaStreamError) -> PyErr {
    match error {
        pdbiox::analysis::SasaStreamError::Surface(pdbiox::surface::SasaError::Spatial(
            pdbiox::SpatialError::Memory(_),
        ))
        | pdbiox::analysis::SasaStreamError::Trajectory(
            pdbiox::traj::TrajectoryError::MemoryLimit { .. },
        ) => PyMemoryError::new_err(error.to_string()),
        _ => crate::core::errors::SasaStreamError::new_err(error.to_string()),
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(visit_shrake_rupley, module)?)?;
    module.add_function(wrap_pyfunction!(collect_shrake_rupley, module)?)?;
    module.add_function(wrap_pyfunction!(sasa_stream, module)?)
}
