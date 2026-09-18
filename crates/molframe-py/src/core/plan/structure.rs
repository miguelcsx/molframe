//! Typed Python views over facade structure operations.

#[path = "structure/config.rs"]
mod config;
#[path = "structure/config_fields.rs"]
mod config_fields;
#[path = "structure/interactions.rs"]
mod interactions;
#[path = "structure/structure_ops.rs"]
mod structure_ops;
#[path = "structure/validation.rs"]
mod validation;
#[path = "structure/values.rs"]
mod values;

pub(crate) use interactions::PyBasePairs;

use super::{Operation, PyPlan, execute_native};
use crate::contract::PyAnalysis;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

pub(super) fn policy_or_default(policy: Option<PyAnalysisPolicy>) -> molframe::AnalysisPolicy {
    policy.map_or_else(molframe::AnalysisPolicy::default, |value| value.inner)
}

pub(super) fn run(
    py: Python<'_>,
    request: &molframe::StructureRequest,
    structure: &PyStructure,
) -> PyResult<PyAnalysis> {
    let plan = PyPlan {
        operations: vec![("result".to_owned(), Operation::Structure(request.clone()))],
    };
    let result = execute_native(&plan, py, Some(structure), None)?;
    let Some(entry) = result.entries.into_iter().next() else {
        return Err(PyValueError::new_err(
            "native structure plan returned no result",
        ));
    };
    match entry.value {
        molframe::PlanValue::Structure(value) => {
            values::to_python(py, *value, Some(structure.structure()))
        }
        _ => Err(PyValueError::new_err(
            "native structure plan returned an incompatible result",
        )),
    }
}

pub(super) fn explain(py: Python<'_>, request: &molframe::StructureRequest) -> PyResult<Py<PyAny>> {
    config::explain(py, request)
}

pub(super) fn to_dict(py: Python<'_>, request: &molframe::StructureRequest) -> PyResult<Py<PyAny>> {
    config::to_dict(py, request)
}

pub(super) fn from_dict(
    config: &Bound<'_, PyDict>,
    expected: &str,
) -> PyResult<molframe::StructureRequest> {
    let request = config::from_dict(config)?;
    if config::kind(&request) != expected {
        return Err(PyValueError::new_err(format!(
            "serialized operation is {}, expected {expected}",
            config::kind(&request)
        )));
    }
    Ok(request)
}

pub(crate) fn typed_request(
    value: &Bound<'_, PyAny>,
) -> PyResult<Option<molframe::StructureRequest>> {
    if let Some(request) = interactions::typed_request(value)? {
        return Ok(Some(request));
    }
    if let Some(request) = structure_ops::typed_request(value)? {
        return Ok(Some(request));
    }
    validation::typed_request(value)
}

pub(crate) fn typed_base_pairs(value: &Bound<'_, PyAny>) -> Option<PyBasePairs> {
    interactions::typed_base_pairs(value)
}

pub(crate) fn serialized_base_pairs(
    py: Python<'_>,
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PyBasePairs>> {
    interactions::serialized_base_pairs(py, config)
}

pub(crate) fn serialized_request(
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<molframe::StructureRequest>> {
    let Some(operation) = config.get_item("operation")? else {
        return Ok(None);
    };
    let operation: String = operation.extract()?;
    if config::is_structure_kind(&operation) {
        return config::from_dict(config).map(Some);
    }
    Ok(None)
}

pub(crate) fn explain_request(
    py: Python<'_>,
    request: &molframe::StructureRequest,
) -> PyResult<Py<PyAny>> {
    explain(py, request)
}

pub(crate) fn dict_request(
    py: Python<'_>,
    request: &molframe::StructureRequest,
) -> PyResult<Py<PyAny>> {
    to_dict(py, request)
}

pub(crate) fn value_to_python(
    py: Python<'_>,
    value: molframe::StructureValue,
    structure: Option<&molframe::Structure>,
) -> PyResult<PyAnalysis> {
    values::to_python(py, value, structure)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    interactions::register(module)?;
    structure_ops::register(module)?;
    validation::register(module)?;
    Ok(())
}
