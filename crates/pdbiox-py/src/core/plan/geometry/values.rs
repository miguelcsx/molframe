//! Native geometry result projection.

use super::model::PyGeometryOperation;
use crate::geometry::{PyAxes, PyDistanceMatrix};
use numpy::IntoPyArray;
use numpy::ndarray::Array1;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub(crate) fn execute_operation(
    py: Python<'_>,
    operation: &PyGeometryOperation,
) -> PyResult<Py<PyAny>> {
    let plan = super::super::PyPlan {
        operations: vec![(
            "result".to_owned(),
            super::super::Operation::Geometry(operation.clone_ref(py)),
        )],
    };
    let native = super::super::execute_native(&plan, py, None, None)?;
    let Some(entry) = native.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native geometry plan returned no result",
        ));
    };
    match entry.value {
        pdbiox::PlanValue::Geometry(value) => value_to_python(py, *value),
        _ => Err(PyValueError::new_err(
            "native geometry plan returned an incompatible result",
        )),
    }
}

pub(crate) fn value_to_python(py: Python<'_>, value: pdbiox::GeometryValue) -> PyResult<Py<PyAny>> {
    match value {
        pdbiox::GeometryValue::Centroid(value) | pdbiox::GeometryValue::CentreOfMass(value) => {
            Ok(value.into_pyobject(py)?.unbind().into_any())
        }
        pdbiox::GeometryValue::RadiusOfGyration(value)
        | pdbiox::GeometryValue::Asphericity(value) => {
            Ok(value.into_pyobject(py)?.unbind().into_any())
        }
        pdbiox::GeometryValue::InertiaTensor(value) => {
            Ok(value.into_pyobject(py)?.unbind().into_any())
        }
        pdbiox::GeometryValue::PrincipalAxes(value)
        | pdbiox::GeometryValue::GyrationAxes(value) => match value {
            Some(value) => Ok(Py::new(py, PyAxes::from(value))?.into_any()),
            None => Ok(py.None()),
        },
        pdbiox::GeometryValue::DistanceMatrix(value) => {
            Ok(Py::new(py, PyDistanceMatrix::from_native(value))?.into_any())
        }
        pdbiox::GeometryValue::Rmsf(value) => {
            Ok(Array1::from_vec(value).into_pyarray(py).unbind().into_any())
        }
    }
}
