//! Serialization and dispatch metadata for physical operations.

use super::{
    PyAnalysisDensityMap, PyCoordinationNumbers, PyLeafletOptions, PyLeaflets, PyLinearDensity,
    PyPhysicalOperation, PyPoreProfile, PyRadialDistribution, PyRadialOptions, PySurfaceContacts,
};
use crate::analysis::{
    PyDensityGridSpec, PyLinearDensityOptions, PyPoreOptions, PySurfaceContactOptions,
};
use crate::graph::PySpatialBackend;
use crate::query::{PyAnalysisPolicy, PySelection};
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

pub(crate) fn policy_or_default(policy: Option<PyAnalysisPolicy>) -> pdbiox::AnalysisPolicy {
    policy.map_or_else(pdbiox::AnalysisPolicy::default, |value| value.inner)
}

pub(crate) fn explain_operation<'py>(
    py: Python<'py>,
    operation: &PyPhysicalOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation.tag())?;
    result.set_item("execution", "native")?;
    result.set_item("gil_released", true)?;
    result.set_item("materializes", "typed analysis result")?;
    result.set_item("complexity", complexity(operation))?;
    Ok(result)
}

pub(crate) fn operation_to_dict<'py>(
    py: Python<'py>,
    operation: &PyPhysicalOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation.tag())?;
    match operation {
        PyPhysicalOperation::RadialDistribution {
            left,
            right,
            options,
            periodic,
            policy,
        } => {
            result.set_item("left", selection_object(py, left)?)?;
            result.set_item("right", selection_object(py, right)?)?;
            result.set_item(
                "options",
                Py::new(py, PyRadialOptions::from_native(*options))?,
            )?;
            result.set_item("periodic", periodic)?;
            result.set_item(
                "policy",
                Py::new(py, PyAnalysisPolicy::from(policy.clone()))?,
            )?;
        }
        PyPhysicalOperation::CoordinationNumbers {
            left,
            right,
            minimum_distance,
            maximum_distance,
            backend,
            periodic,
            policy,
        } => {
            result.set_item("left", selection_object(py, left)?)?;
            result.set_item("right", selection_object(py, right)?)?;
            result.set_item("shell", (*minimum_distance, *maximum_distance))?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
            result.set_item("periodic", periodic)?;
            result.set_item(
                "policy",
                Py::new(py, PyAnalysisPolicy::from(policy.clone()))?,
            )?;
        }
        PyPhysicalOperation::Leaflets {
            sites,
            options,
            periodic,
            policy,
        } => {
            result.set_item("sites", selection_object(py, sites)?)?;
            result.set_item(
                "options",
                Py::new(py, PyLeafletOptions::from_native(*options))?,
            )?;
            result.set_item("periodic", periodic)?;
            result.set_item(
                "policy",
                Py::new(py, PyAnalysisPolicy::from(policy.clone()))?,
            )?;
        }
        PyPhysicalOperation::LinearDensity { .. }
        | PyPhysicalOperation::DensityMap { .. }
        | PyPhysicalOperation::PoreProfile { .. }
        | PyPhysicalOperation::SurfaceContacts { .. } => {
            array_operation_to_dict(py, operation, &result)?;
        }
    }
    Ok(result)
}

fn array_operation_to_dict(
    py: Python<'_>,
    operation: &PyPhysicalOperation,
    result: &Bound<'_, PyDict>,
) -> PyResult<()> {
    match operation {
        PyPhysicalOperation::LinearDensity {
            weights,
            options,
            policy,
        } => {
            result.set_item("weights", weights.bind(py))?;
            result.set_item(
                "options",
                Py::new(py, PyLinearDensityOptions::from_native(*options))?,
            )?;
            set_policy(py, result, policy)?;
        }
        PyPhysicalOperation::DensityMap {
            weights,
            spec,
            policy,
        } => {
            result.set_item("weights", weights.bind(py))?;
            result.set_item("spec", Py::new(py, PyDensityGridSpec::from_native(*spec))?)?;
            set_policy(py, result, policy)?;
        }
        PyPhysicalOperation::PoreProfile {
            radii,
            options,
            policy,
        } => {
            result.set_item("radii", radii.bind(py))?;
            result.set_item(
                "options",
                Py::new(py, PyPoreOptions::from_native(*options))?,
            )?;
            set_policy(py, result, policy)?;
        }
        PyPhysicalOperation::SurfaceContacts {
            radii,
            tolerance,
            probe,
            density,
            minimum_area,
            backend,
            policy,
        } => {
            result.set_item("radii", radii.bind(py))?;
            result.set_item(
                "options",
                Py::new(
                    py,
                    PySurfaceContactOptions::from_native_parts(
                        *tolerance,
                        *probe,
                        *density,
                        *minimum_area,
                        *backend,
                    ),
                )?,
            )?;
            set_policy(py, result, policy)?;
        }
        _ => return Err(PyValueError::new_err("not an array physical operation")),
    }
    Ok(())
}

fn set_policy(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    policy: &pdbiox::AnalysisPolicy,
) -> PyResult<()> {
    result.set_item(
        "policy",
        Py::new(py, PyAnalysisPolicy::from(policy.clone()))?,
    )
}

pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<PyPhysicalOperation> {
    let operation: String = required(config, "operation")?.extract()?;
    let policy = policy_from_dict(config)?;
    match operation.as_str() {
        "radial_distribution" => Ok(PyPhysicalOperation::RadialDistribution {
            left: required(config, "left")?.extract::<PySelection>()?.inner,
            right: required(config, "right")?.extract::<PySelection>()?.inner,
            options: required(config, "options")?
                .extract::<PyRadialOptions>()?
                .native(),
            periodic: optional(config, "periodic")?.map_or(Ok(false), |value| value.extract())?,
            policy,
        }),
        "coordination_numbers" => {
            let shell: (f32, f32) = required(config, "shell")?.extract()?;
            Ok(PyPhysicalOperation::CoordinationNumbers {
                left: required(config, "left")?.extract::<PySelection>()?.inner,
                right: required(config, "right")?.extract::<PySelection>()?.inner,
                minimum_distance: shell.0,
                maximum_distance: shell.1,
                backend: required(config, "backend")?
                    .extract::<PySpatialBackend>()?
                    .into(),
                periodic: optional(config, "periodic")?
                    .map_or(Ok(false), |value| value.extract())?,
                policy,
            })
        }
        "leaflets" => Ok(PyPhysicalOperation::Leaflets {
            sites: required(config, "sites")?.extract::<PySelection>()?.inner,
            options: required(config, "options")?
                .extract::<PyLeafletOptions>()?
                .native(),
            periodic: optional(config, "periodic")?.map_or(Ok(false), |value| value.extract())?,
            policy,
        }),
        "linear_density" => Ok(PyPhysicalOperation::LinearDensity {
            weights: required(config, "weights")?.unbind(),
            options: required(config, "options")?
                .extract::<PyLinearDensityOptions>()?
                .native(),
            policy,
        }),
        "density_map" => Ok(PyPhysicalOperation::DensityMap {
            weights: required(config, "weights")?.unbind(),
            spec: required(config, "spec")?
                .extract::<PyDensityGridSpec>()?
                .native(),
            policy,
        }),
        "pore_profile" => Ok(PyPhysicalOperation::PoreProfile {
            radii: required(config, "radii")?.unbind(),
            options: required(config, "options")?
                .extract::<PyPoreOptions>()?
                .native(),
            policy,
        }),
        "surface_contacts" => {
            let options = required(config, "options")?.extract::<PySurfaceContactOptions>()?;
            let (tolerance, probe, density, minimum_area, backend) = options.native_parts();
            Ok(PyPhysicalOperation::SurfaceContacts {
                radii: required(config, "radii")?.unbind(),
                tolerance,
                probe,
                density,
                minimum_area,
                backend,
                policy,
            })
        }
        _ => Err(PyValueError::new_err(format!(
            "unknown physical operation {operation:?}"
        ))),
    }
}

pub(crate) fn typed_request(value: &Bound<'_, PyAny>) -> Option<PyPhysicalOperation> {
    if let Ok(value) = value.extract::<PyRef<'_, PyRadialDistribution>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyCoordinationNumbers>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyLeaflets>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyLinearDensity>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyAnalysisDensityMap>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyPoreProfile>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PySurfaceContacts>>() {
        return Some(value.operation.clone_ref(value.py()));
    }
    None
}

pub(crate) fn serialized_request(
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PyPhysicalOperation>> {
    let Some(operation) = config.get_item("operation")? else {
        return Ok(None);
    };
    let operation: &str = operation.extract()?;
    if matches!(
        operation,
        "radial_distribution"
            | "coordination_numbers"
            | "leaflets"
            | "linear_density"
            | "density_map"
            | "pore_profile"
            | "surface_contacts"
    ) {
        return from_dict(config).map(Some);
    }
    Ok(None)
}

fn policy_from_dict(config: &Bound<'_, PyDict>) -> PyResult<pdbiox::AnalysisPolicy> {
    match config.get_item("policy")? {
        Some(value) if !value.is_none() => Ok(value.extract::<PyAnalysisPolicy>()?.inner),
        _ => Ok(pdbiox::AnalysisPolicy::default()),
    }
}

fn complexity(operation: &PyPhysicalOperation) -> &'static str {
    match operation {
        PyPhysicalOperation::RadialDistribution { .. }
        | PyPhysicalOperation::CoordinationNumbers { .. } => "O(atoms + local_neighbours)",
        PyPhysicalOperation::Leaflets { .. } => "O(sites + local_neighbours + sites·α(sites))",
        PyPhysicalOperation::LinearDensity { .. } => "O(atoms + bins)",
        PyPhysicalOperation::DensityMap { .. } => "O(atoms + grid_cells)",
        PyPhysicalOperation::PoreProfile { .. } => "O(atoms + samples × local_neighbours)",
        PyPhysicalOperation::SurfaceContacts { .. } => {
            "O(atoms + local_neighbours × surface_samples)"
        }
    }
}

fn selection_object(py: Python<'_>, selection: &pdbiox::AtomSelection) -> PyResult<Py<PyAny>> {
    Ok(Py::new(
        py,
        PySelection {
            inner: selection.clone(),
        },
    )?
    .into_any())
}

fn required<'py>(config: &'py Bound<'py, PyDict>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    config
        .get_item(name)?
        .ok_or_else(|| PyKeyError::new_err(format!("missing operation field: {name}")))
}

fn optional<'py>(
    config: &'py Bound<'py, PyDict>,
    name: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    config.get_item(name)
}

pub(super) fn expected(operation: &str) -> PyErr {
    PyValueError::new_err(format!("serialized operation is not {operation}"))
}
