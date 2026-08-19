//! Native bindings for formatted AMBER restart and ASCII trajectories.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, AmberError, pyo3::exceptions::PyValueError);

#[pyclass(name = "AmberRestartLayout", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAmberRestartLayout {
    Auto,
    Coordinates,
    CoordinatesBox3,
    CoordinatesBox6,
    CoordinatesVelocities,
    CoordinatesVelocitiesBox3,
    CoordinatesVelocitiesBox6,
}

impl From<PyAmberRestartLayout> for pdbiox::traj::AmberRestartLayout {
    fn from(value: PyAmberRestartLayout) -> Self {
        match value {
            PyAmberRestartLayout::Auto => Self::Auto,
            PyAmberRestartLayout::Coordinates => Self::Coordinates,
            PyAmberRestartLayout::CoordinatesBox3 => Self::CoordinatesBox3,
            PyAmberRestartLayout::CoordinatesBox6 => Self::CoordinatesBox6,
            PyAmberRestartLayout::CoordinatesVelocities => Self::CoordinatesVelocities,
            PyAmberRestartLayout::CoordinatesVelocitiesBox3 => Self::CoordinatesVelocitiesBox3,
            PyAmberRestartLayout::CoordinatesVelocitiesBox6 => Self::CoordinatesVelocitiesBox6,
        }
    }
}

impl From<pdbiox::traj::AmberRestartLayout> for PyAmberRestartLayout {
    fn from(value: pdbiox::traj::AmberRestartLayout) -> Self {
        match value {
            pdbiox::traj::AmberRestartLayout::Auto => Self::Auto,
            pdbiox::traj::AmberRestartLayout::Coordinates => Self::Coordinates,
            pdbiox::traj::AmberRestartLayout::CoordinatesBox3 => Self::CoordinatesBox3,
            pdbiox::traj::AmberRestartLayout::CoordinatesBox6 => Self::CoordinatesBox6,
            pdbiox::traj::AmberRestartLayout::CoordinatesVelocities => Self::CoordinatesVelocities,
            pdbiox::traj::AmberRestartLayout::CoordinatesVelocitiesBox3 => {
                Self::CoordinatesVelocitiesBox3
            }
            pdbiox::traj::AmberRestartLayout::CoordinatesVelocitiesBox6 => {
                Self::CoordinatesVelocitiesBox6
            }
        }
    }
}

#[pyclass(name = "AmberRestart", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberRestart {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) layout: PyAmberRestartLayout,
    #[pyo3(get)]
    pub(crate) timestep: PyTimestep,
}

#[pymethods]
impl PyAmberRestart {
    #[new]
    fn new(title: String, layout: PyAmberRestartLayout, timestep: PyTimestep) -> Self {
        Self {
            title,
            layout,
            timestep,
        }
    }
}

impl TryFrom<PyAmberRestart> for pdbiox::traj::AmberRestart {
    type Error = PyErr;

    fn try_from(value: PyAmberRestart) -> Result<Self, Self::Error> {
        Ok(Self {
            title: value.title.into(),
            layout: value.layout.into(),
            timestep: value.timestep.try_into()?,
        })
    }
}

impl TryFrom<pdbiox::traj::AmberRestart> for PyAmberRestart {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::AmberRestart) -> Result<Self, Self::Error> {
        Ok(Self {
            title: value.title.into(),
            layout: value.layout.into(),
            timestep: value.timestep.try_into()?,
        })
    }
}

#[pyfunction]
#[pyo3(signature = (text, layout))]
pub(crate) fn parse_amber_restart(
    text: &str,
    layout: PyAmberRestartLayout,
) -> PyResult<PyTimestep> {
    pdbiox::traj::parse_amber_restart(text, layout.into())
        .map_err(|error| AmberError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
#[pyo3(signature = (text, layout))]
pub(crate) fn parse_amber_restart_record(
    text: &str,
    layout: PyAmberRestartLayout,
) -> PyResult<PyAmberRestart> {
    pdbiox::traj::parse_amber_restart_record(text, layout.into())
        .map_err(|error| AmberError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_amber_restart(record: PyAmberRestart) -> PyResult<String> {
    let record = record.try_into()?;
    pdbiox::traj::write_amber_restart(&record)
        .map_err(|error| AmberError::new_err(error.to_string()))
}

#[pyfunction]
#[pyo3(signature = (text, atom_count, periodic_box))]
pub(crate) fn parse_amber_ascii_trajectory(
    text: &str,
    atom_count: usize,
    periodic_box: bool,
) -> PyResult<Vec<PyTimestep>> {
    pdbiox::traj::parse_amber_ascii_trajectory(text, atom_count, periodic_box)
        .map_err(|error| AmberError::new_err(error.to_string()))?
        .into_iter()
        .map(TryInto::try_into)
        .collect()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("AmberError", module.py().get_type::<AmberError>())?;
    module.add_class::<PyAmberRestartLayout>()?;
    module.add_class::<PyAmberRestart>()?;
    module.add_function(wrap_pyfunction!(parse_amber_restart, module)?)?;
    module.add_function(wrap_pyfunction!(parse_amber_restart_record, module)?)?;
    module.add_function(wrap_pyfunction!(write_amber_restart, module)?)?;
    module.add_function(wrap_pyfunction!(parse_amber_ascii_trajectory, module)?)?;
    Ok(())
}
