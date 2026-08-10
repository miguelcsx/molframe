//! One-call Python entry points into the Rust trajectory dispatcher.

use super::{PyTrajectory, PyTrajectoryFormat, PyTrajectoryUnits, PyTrajectoryWriteOptions};
use pdbiox::traj::{
    GsdOptions, H5mdOptions, TrajectoryFormat, TrajectoryReadOptions, TrajectoryWriteOptions,
    read_trajectory, write_trajectory,
};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

const UNKNOWN_OUTPUT_FORMAT: &str = "trajectory output format cannot be inferred from the path";

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
        let path = path.into_boxed_path();
        let units = units.map(|units| units.inner.clone());
        let mut options = TrajectoryReadOptions {
            format: format.map(Into::into),
            gsd: length_to_angstrom
                .map(GsdOptions::new)
                .transpose()
                .map_err(value_error)?,
            ..TrajectoryReadOptions::default()
        };
        apply_h5md_options(&mut options.h5md, particle_group, units);
        read_trajectory(&path, &options)
            .map_err(io_error)
            .and_then(PyTrajectory::from_data)
    }

    #[pyo3(signature = (path, *, options=None))]
    fn write(
        &self,
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
        write_trajectory(&path, &self.to_data(format), &writer).map_err(io_error)
    }
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
