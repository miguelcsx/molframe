//! Python projections for the reusable trajectory reader data model.

use crate::crystallography::PyUnitCell;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "FrameValue", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFrameValue {
    #[pyo3(get)]
    pub(crate) kind: String,
    #[pyo3(get)]
    pub(crate) float_value: Option<f64>,
    #[pyo3(get)]
    pub(crate) integer_value: Option<i64>,
    #[pyo3(get)]
    pub(crate) text_value: Option<String>,
    #[pyo3(get)]
    pub(crate) floats_value: Option<Vec<f64>>,
}

#[pymethods]
impl PyFrameValue {
    #[staticmethod]
    fn float(value: f64) -> Self {
        Self {
            kind: "float".into(),
            float_value: Some(value),
            integer_value: None,
            text_value: None,
            floats_value: None,
        }
    }

    #[staticmethod]
    fn integer(value: i64) -> Self {
        Self {
            kind: "integer".into(),
            float_value: None,
            integer_value: Some(value),
            text_value: None,
            floats_value: None,
        }
    }

    #[staticmethod]
    fn text(value: String) -> Self {
        Self {
            kind: "text".into(),
            float_value: None,
            integer_value: None,
            text_value: Some(value),
            floats_value: None,
        }
    }

    #[staticmethod]
    fn floats(value: Vec<f64>) -> Self {
        Self {
            kind: "floats".into(),
            float_value: None,
            integer_value: None,
            text_value: None,
            floats_value: Some(value),
        }
    }
}

impl TryFrom<molframe::traj::FrameValue> for PyFrameValue {
    type Error = PyErr;

    fn try_from(value: molframe::traj::FrameValue) -> Result<Self, Self::Error> {
        match value {
            molframe::traj::FrameValue::Float(value) => Ok(Self::float(value)),
            molframe::traj::FrameValue::Integer(value) => Ok(Self::integer(value)),
            molframe::traj::FrameValue::Text(value) => Ok(Self::text(value.into())),
            molframe::traj::FrameValue::Floats(value) => Ok(Self::floats(value)),
            _ => Err(pyo3::exceptions::PyValueError::new_err(
                "native frame-value variant is not supported by this Python version",
            )),
        }
    }
}

impl TryFrom<PyFrameValue> for molframe::traj::FrameValue {
    type Error = PyErr;

    fn try_from(value: PyFrameValue) -> Result<Self, Self::Error> {
        match value.kind.as_str() {
            "float" => value
                .float_value
                .map(Self::Float)
                .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("float value is absent")),
            "integer" => value
                .integer_value
                .map(Self::Integer)
                .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("integer value is absent")),
            "text" => value
                .text_value
                .map(|value| Self::Text(value.into()))
                .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("text value is absent")),
            "floats" => value
                .floats_value
                .map(Self::Floats)
                .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("float vector is absent")),
            _ => Err(pyo3::exceptions::PyValueError::new_err(
                "unknown frame-value kind",
            )),
        }
    }
}

#[pyclass(name = "Timestep", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTimestep {
    #[pyo3(get)]
    pub(crate) frame: usize,
    #[pyo3(get)]
    pub(crate) time: Option<f64>,
    #[pyo3(get)]
    pub(crate) dt: Option<f64>,
    #[pyo3(get)]
    pub(crate) positions: Vec<[f32; 3]>,
    #[pyo3(get)]
    pub(crate) velocities: Option<Vec<[f32; 3]>>,
    #[pyo3(get)]
    pub(crate) forces: Option<Vec<[f32; 3]>>,
    #[pyo3(get)]
    pub(crate) cell: Option<PyUnitCell>,
    #[pyo3(get)]
    pub(crate) data: Vec<(String, PyFrameValue)>,
}

#[pymethods]
impl PyTimestep {
    #[new]
    #[pyo3(signature = (positions, *, frame=0, time=None, dt=None, velocities=None, forces=None, cell=None, data=None))]
    fn new(
        positions: Vec<[f32; 3]>,
        frame: usize,
        time: Option<f64>,
        dt: Option<f64>,
        velocities: Option<Vec<[f32; 3]>>,
        forces: Option<Vec<[f32; 3]>>,
        cell: Option<PyRef<'_, PyUnitCell>>,
        data: Option<Vec<(String, PyFrameValue)>>,
    ) -> PyResult<Self> {
        if velocities
            .as_ref()
            .is_some_and(|values| values.len() != positions.len())
            || forces
                .as_ref()
                .is_some_and(|values| values.len() != positions.len())
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "timestep streams must match the position atom count",
            ));
        }
        Ok(Self {
            frame,
            time,
            dt,
            positions,
            velocities,
            forces,
            cell: cell.map(|cell| cell.clone()),
            data: match data {
                Some(data) => data,
                None => Vec::new(),
            },
        })
    }
}

impl TryFrom<PyTimestep> for molframe::traj::Timestep {
    type Error = PyErr;

    fn try_from(value: PyTimestep) -> Result<Self, Self::Error> {
        let data = value
            .data
            .into_iter()
            .map(|(key, value)| value.try_into().map(|value| (key.into(), value)))
            .collect::<Result<BTreeMap<_, _>, PyErr>>()?;
        Ok(Self {
            frame: value.frame,
            time: value.time,
            dt: value.dt,
            positions: value.positions,
            velocities: value.velocities,
            forces: value.forces,
            cell: value.cell.map(|cell| cell.cell),
            data,
        })
    }
}

impl TryFrom<molframe::traj::Timestep> for PyTimestep {
    type Error = PyErr;

    fn try_from(value: molframe::traj::Timestep) -> Result<Self, Self::Error> {
        let cell = value.cell.map(PyUnitCell::from_native).transpose()?;
        let data = value
            .data
            .into_iter()
            .map(|(key, value)| value.try_into().map(|value| (key.into(), value)))
            .collect::<Result<_, PyErr>>()?;
        Ok(Self {
            frame: value.frame,
            time: value.time,
            dt: value.dt,
            positions: value.positions,
            velocities: value.velocities,
            forces: value.forces,
            cell,
            data,
        })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFrameValue>()?;
    module.add_class::<PyTimestep>()?;
    Ok(())
}
