//! Typed and serialized operation-node dispatch for the Python plan.

use super::coordinates::OperationKind;
use super::{
    CoordinateMetric, Operation, PyAny, PyComparison, PyContacts, PyGdtHa, PyGdtTs, PyLddt, PyRef,
    PyRmsd, PyTmScore, chemistry, geometry, physical, required, selection, spatial, structure,
    surface, trajectory,
};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

pub(super) fn operation_from_value(value: Bound<'_, PyAny>, name: &str) -> PyResult<Operation> {
    if let Some(operation) = typed_operation(&value)? {
        return Ok(operation);
    }
    if let Ok(config) = value.cast::<PyDict>()
        && let Some(operation) = serialized_operation(value.py(), config)?
    {
        return Ok(operation);
    }
    Err(PyTypeError::new_err(format!(
        "operation {name:?} must be a typed native operation or a serialized operation"
    )))
}

/// Requires a serialized configuration to belong to one concrete operation.
///
/// Public `from_dict` constructors are intentionally strict: accepting another
/// operation's fields would create a valid-looking object with different
/// semantics from the configuration that produced it.
pub(super) fn require_operation_tag(config: &Bound<'_, PyDict>, expected: &str) -> PyResult<()> {
    let actual: String = required(config, "operation")?.extract()?;
    if actual == expected {
        Ok(())
    } else {
        Err(PyValueError::new_err(format!(
            "expected serialized operation {expected:?}, found {actual:?}"
        )))
    }
}

fn typed_operation(value: &Bound<'_, PyAny>) -> PyResult<Option<Operation>> {
    if let Ok(value) = value.extract::<PyRef<'_, PyContacts>>() {
        return Ok(Some(Operation::Contacts(value.clone())));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyRmsd>>() {
        return Ok(Some(Operation::Rmsd(PyRmsd {
            mobile: value.mobile.clone_ref(value.py()),
            reference: value.reference.clone_ref(value.py()),
        })));
    }
    if let Some(operation) = comparison_operation(value) {
        return Ok(Some(Operation::Comparison(operation)));
    }
    if let Ok(value) = value.extract::<PyRef<'_, chemistry::PyInferBonds>>() {
        return Ok(Some(Operation::BondInference(value.clone())));
    }
    if let Some(operation) = selection::typed_request(value) {
        return Ok(Some(Operation::Selection(operation)));
    }
    if let Some(operation) = geometry::typed_request(value)? {
        return Ok(Some(Operation::Geometry(operation)));
    }
    if let Some(operation) = trajectory::typed_request(value)? {
        return Ok(Some(Operation::Trajectory(operation)));
    }
    if let Some(operation) = spatial::typed_request(value) {
        return Ok(Some(Operation::Spatial(operation)));
    }
    if let Some(operation) = physical::typed_request(value) {
        return Ok(Some(Operation::Physical(operation)));
    }
    if let Some(operation) = surface::typed_request(value) {
        return Ok(Some(Operation::Surface(operation)));
    }
    if let Some(operation) = structure::typed_base_pairs(value) {
        return Ok(Some(Operation::BasePairs(operation)));
    }
    if let Some(request) = structure::typed_request(value)? {
        return Ok(Some(Operation::Structure(request)));
    }
    Ok(None)
}

fn comparison_operation(value: &Bound<'_, PyAny>) -> Option<PyComparison> {
    if let Ok(value) = value.extract::<PyRef<'_, PyLddt>>() {
        return Some(PyComparison {
            mobile: value.mobile.clone_ref(value.py()),
            reference: value.reference.clone_ref(value.py()),
            metric: CoordinateMetric::Lddt {
                inclusion_radius: value.inclusion_radius,
            },
        });
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyTmScore>>() {
        return Some(PyComparison {
            mobile: value.mobile.clone_ref(value.py()),
            reference: value.reference.clone_ref(value.py()),
            metric: CoordinateMetric::TmScore,
        });
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyGdtTs>>() {
        return Some(PyComparison {
            mobile: value.mobile.clone_ref(value.py()),
            reference: value.reference.clone_ref(value.py()),
            metric: CoordinateMetric::GdtTs,
        });
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyGdtHa>>() {
        return Some(PyComparison {
            mobile: value.mobile.clone_ref(value.py()),
            reference: value.reference.clone_ref(value.py()),
            metric: CoordinateMetric::GdtHa,
        });
    }
    None
}

fn serialized_operation(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Option<Operation>> {
    if let Some(operation) = selection::serialized_request(py, config)? {
        return Ok(Some(Operation::Selection(operation)));
    }
    if let Some(operation) = geometry::serialized_request(config)? {
        return Ok(Some(Operation::Geometry(operation)));
    }
    if let Some(operation) = trajectory::serialized_request(config)? {
        return Ok(Some(Operation::Trajectory(operation)));
    }
    if let Some(operation) = spatial::serialized_request(py, config)? {
        return Ok(Some(Operation::Spatial(operation)));
    }
    if let Some(operation) = physical::serialized_request(config)? {
        return Ok(Some(Operation::Physical(operation)));
    }
    if let Some(operation) = surface::serialized_request(py, config)? {
        return Ok(Some(Operation::Surface(operation)));
    }
    if let Some(operation) = chemistry::serialized_request(config)? {
        return Ok(Some(Operation::BondInference(operation)));
    }
    if let Some(operation) = structure::serialized_base_pairs(py, config)? {
        return Ok(Some(Operation::BasePairs(operation)));
    }
    if let Some(request) = structure::serialized_request(config)? {
        return Ok(Some(Operation::Structure(request)));
    }
    let operation: String = required(config, "operation")?.extract()?;
    let result = match OperationKind::parse(&operation) {
        Some(OperationKind::Contacts) => Operation::Contacts(PyContacts::from_dict(config)?),
        Some(OperationKind::Rmsd) => Operation::Rmsd(PyRmsd::from_dict(config)?),
        Some(OperationKind::Lddt) => comparison_from_config(config, OperationKind::Lddt)?,
        Some(OperationKind::TmScore) => comparison_from_config(config, OperationKind::TmScore)?,
        Some(OperationKind::GdtTs) => comparison_from_config(config, OperationKind::GdtTs)?,
        Some(OperationKind::GdtHa) => comparison_from_config(config, OperationKind::GdtHa)?,
        None => {
            return Err(PyValueError::new_err(format!(
                "unknown serialized operation {operation:?}"
            )));
        }
    };
    Ok(Some(result))
}

fn comparison_from_config(config: &Bound<'_, PyDict>, kind: OperationKind) -> PyResult<Operation> {
    let (mobile, reference, metric) = match kind {
        OperationKind::Lddt => {
            let value = PyLddt::from_dict(config)?;
            (
                value.mobile,
                value.reference,
                CoordinateMetric::Lddt {
                    inclusion_radius: value.inclusion_radius,
                },
            )
        }
        OperationKind::TmScore => {
            let value = PyTmScore::from_dict(config)?;
            (value.mobile, value.reference, CoordinateMetric::TmScore)
        }
        OperationKind::GdtTs => {
            let value = PyGdtTs::from_dict(config)?;
            (value.mobile, value.reference, CoordinateMetric::GdtTs)
        }
        OperationKind::GdtHa => {
            let value = PyGdtHa::from_dict(config)?;
            (value.mobile, value.reference, CoordinateMetric::GdtHa)
        }
        OperationKind::Contacts | OperationKind::Rmsd => {
            return Err(PyValueError::new_err("not a comparison operation"));
        }
    };
    Ok(Operation::Comparison(PyComparison {
        mobile,
        reference,
        metric,
    }))
}

#[cfg(test)]
#[path = "dispatch_tests.rs"]
mod tests;
