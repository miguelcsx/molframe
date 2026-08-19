//! One-call Python entry points into the Rust trajectory dispatcher.

use super::{
    PyAmberRestartLayout, PyTrajectory, PyTrajectoryFormat, PyTrajectoryUnits,
    PyTrajectoryWriteOptions,
};
use pdbiox::traj::{
    AmberAsciiReadOptions, GsdOptions, H5mdOptions, TrajectoryFormat, TrajectoryReadOptions,
    TrajectoryWriteOptions, read_trajectory as native_read_trajectory,
    write_trajectory as native_write_trajectory,
};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

const UNKNOWN_OUTPUT_FORMAT: &str = "trajectory output format cannot be inferred from the path";

#[pyclass(name = "AmberAsciiReadOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAmberAsciiReadOptions {
    #[pyo3(get)]
    pub(crate) atom_count: usize,
    #[pyo3(get)]
    pub(crate) periodic_box: bool,
}

#[pymethods]
impl PyAmberAsciiReadOptions {
    #[new]
    #[pyo3(signature = (atom_count, *, periodic_box=false))]
    fn new(atom_count: usize, periodic_box: bool) -> PyResult<Self> {
        if atom_count == 0 {
            return Err(PyValueError::new_err(
                "atom_count must be greater than zero",
            ));
        }
        Ok(Self {
            atom_count,
            periodic_box,
        })
    }
}

#[pyclass(name = "TrajectoryReadOptions", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTrajectoryReadOptions {
    pub(crate) inner: TrajectoryReadOptions,
}

#[pymethods]
impl PyTrajectoryReadOptions {
    #[new]
    #[pyo3(signature = (*, format=None, length_to_angstrom=None, particle_group=None, units=None, amber_restart_layout=None, amber_ascii=None))]
    fn new(
        format: Option<PyTrajectoryFormat>,
        length_to_angstrom: Option<f64>,
        particle_group: Option<String>,
        units: Option<PyRef<'_, PyTrajectoryUnits>>,
        amber_restart_layout: Option<PyAmberRestartLayout>,
        amber_ascii: Option<PyRef<'_, PyAmberAsciiReadOptions>>,
    ) -> PyResult<Self> {
        let mut inner = TrajectoryReadOptions {
            format: format.map(Into::into),
            gsd: length_to_angstrom
                .map(GsdOptions::new)
                .transpose()
                .map_err(value_error)?,
            amber_restart_layout: match amber_restart_layout {
                Some(layout) => layout.into(),
                None => pdbiox::traj::AmberRestartLayout::Auto,
            },
            amber_ascii: amber_ascii.map(|options| AmberAsciiReadOptions {
                atom_count: options.atom_count,
                periodic_box: options.periodic_box,
            }),
            ..TrajectoryReadOptions::default()
        };
        apply_h5md_options(
            &mut inner.h5md,
            particle_group,
            units.map(|value| value.inner.clone()),
        );
        Ok(Self { inner })
    }

    #[staticmethod]
    fn defaults() -> Self {
        Self {
            inner: TrajectoryReadOptions::default(),
        }
    }
}

#[pymethods]
impl PyTrajectory {
    #[staticmethod]
    #[pyo3(signature = (path, *, format=None, length_to_angstrom=None, particle_group=None, units=None))]
    fn read(
        path: PathBuf,
        format: Option<PyTrajectoryFormat>,
        length_to_angstrom: Option<f64>,
        particle_group: Option<String>,
        units: Option<PyRef<'_, PyTrajectoryUnits>>,
    ) -> PyResult<Self> {
        let options = PyTrajectoryReadOptions::new(
            format,
            length_to_angstrom,
            particle_group,
            units,
            None,
            None,
        )?;
        read_native(path, &options.inner)
    }

    #[pyo3(signature = (path, *, options=None))]
    fn write(
        &self,
        py: Python<'_>,
        path: PathBuf,
        options: Option<PyRef<'_, PyTrajectoryWriteOptions>>,
    ) -> PyResult<()> {
        let path = path.into_boxed_path();
        let writer = options.map_or_else(TrajectoryWriteOptions::default, |options| {
            options.inner.clone()
        });
        let format = writer
            .format
            .or_else(|| TrajectoryFormat::infer(&path))
            .or(self.format)
            .ok_or_else(|| PyValueError::new_err(UNKNOWN_OUTPUT_FORMAT))?;
        native_write_trajectory(&path, &self.to_data(py, format)?, &writer).map_err(io_error)
    }
}

#[pyfunction]
#[pyo3(signature = (path, *, options=None))]
pub(crate) fn read_trajectory(
    path: PathBuf,
    options: Option<PyRef<'_, PyTrajectoryReadOptions>>,
) -> PyResult<PyTrajectory> {
    let options = options.map_or_else(TrajectoryReadOptions::default, |value| value.inner.clone());
    read_native(path, &options)
}

#[pyfunction]
#[pyo3(signature = (path, trajectory, *, options=None))]
pub(crate) fn write_trajectory(
    py: Python<'_>,
    path: PathBuf,
    trajectory: &PyTrajectory,
    options: Option<PyRef<'_, PyTrajectoryWriteOptions>>,
) -> PyResult<()> {
    let writer = options.map_or_else(TrajectoryWriteOptions::default, |value| value.inner.clone());
    let format = writer
        .format
        .or_else(|| TrajectoryFormat::infer(&path))
        .or(trajectory.format)
        .ok_or_else(|| PyValueError::new_err(UNKNOWN_OUTPUT_FORMAT))?;
    native_write_trajectory(&path, &trajectory.to_data(py, format)?, &writer).map_err(io_error)
}

fn read_native(path: PathBuf, options: &TrajectoryReadOptions) -> PyResult<PyTrajectory> {
    let path = path.into_boxed_path();
    native_read_trajectory(&path, options)
        .map_err(io_error)
        .and_then(PyTrajectory::from_data)
}

fn apply_h5md_options(
    options: &mut H5mdOptions,
    particle_group: Option<String>,
    units: Option<pdbiox::traj::H5mdUnitSystem>,
) {
    if let Some(particle_group) = particle_group {
        options.particle_group = particle_group;
    }
    if let Some(units) = units {
        options.units = units;
    }
}

fn io_error(error: pdbiox::traj::TrajectoryIoError) -> PyErr {
    match error {
        pdbiox::traj::TrajectoryIoError::Io(error) => PyOSError::new_err(error.to_string()),
        error => PyValueError::new_err(error.to_string()),
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
