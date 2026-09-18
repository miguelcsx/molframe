//! Serialization and execution metadata for structure operations.

use crate::analysis::{
    PyBasePairOptions, PyCationPiOptions, PyDsspOptions, PyGnmOptions, PyPiStackingOptions,
    PyPlanarityOptions, PyWaterBridgeOptions,
};
use crate::chemistry::PyRadiusSet;
use crate::graph::PySpatialBackend;
use crate::query::{PyAnalysisPolicy, PySelection};
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::sync::Arc;

pub(super) fn kind(request: &molframe::StructureRequest) -> &'static str {
    match request {
        molframe::StructureRequest::BasePairs { .. } => "base_pairs",
        molframe::StructureRequest::HydrogenBonds { .. } => "hydrogen_bonds",
        molframe::StructureRequest::SaltBridges { .. } => "salt_bridges",
        molframe::StructureRequest::PiStacking { .. } => "pi_stacking",
        molframe::StructureRequest::CationPi { .. } => "cation_pi",
        molframe::StructureRequest::WaterBridges { .. } => "water_bridges",
        molframe::StructureRequest::ContactMap { .. } => "contact_map",
        molframe::StructureRequest::ChainInterface { .. } => "chain_interface",
        molframe::StructureRequest::SecondaryStructure { .. } => "secondary_structure",
        molframe::StructureRequest::HalfSphereExposure { .. } => "half_sphere_exposure",
        molframe::StructureRequest::NucleicTorsions { .. } => "nucleic_torsions",
        molframe::StructureRequest::GaussianNetworkModel { .. } => "gaussian_network_model",
        molframe::StructureRequest::Clashes { .. } => "clashes",
        molframe::StructureRequest::BondLengthDeviations { .. } => "bond_length_deviations",
        molframe::StructureRequest::CisPeptides { .. } => "cis_peptides",
        molframe::StructureRequest::Planarity { .. } => "planarity",
        molframe::StructureRequest::Quality { .. } => "quality",
        molframe::StructureRequest::Valence { .. } => "valence",
        molframe::StructureRequest::Completeness { .. } => "completeness",
    }
}

pub(super) fn is_structure_kind(value: &str) -> bool {
    matches!(
        value,
        "base_pairs"
            | "hydrogen_bonds"
            | "salt_bridges"
            | "pi_stacking"
            | "cation_pi"
            | "water_bridges"
            | "contact_map"
            | "chain_interface"
            | "secondary_structure"
            | "half_sphere_exposure"
            | "nucleic_torsions"
            | "gaussian_network_model"
            | "clashes"
            | "bond_length_deviations"
            | "cis_peptides"
            | "planarity"
            | "quality"
            | "valence"
            | "completeness"
    )
}

fn policy(request: &molframe::StructureRequest) -> &molframe::AnalysisPolicy {
    match request {
        molframe::StructureRequest::BasePairs { policy, .. }
        | molframe::StructureRequest::HydrogenBonds { policy, .. }
        | molframe::StructureRequest::SaltBridges { policy, .. }
        | molframe::StructureRequest::PiStacking { policy, .. }
        | molframe::StructureRequest::CationPi { policy, .. }
        | molframe::StructureRequest::WaterBridges { policy, .. }
        | molframe::StructureRequest::ContactMap { policy, .. }
        | molframe::StructureRequest::ChainInterface { policy, .. }
        | molframe::StructureRequest::SecondaryStructure { policy, .. }
        | molframe::StructureRequest::HalfSphereExposure { policy, .. }
        | molframe::StructureRequest::NucleicTorsions { policy }
        | molframe::StructureRequest::GaussianNetworkModel { policy, .. }
        | molframe::StructureRequest::Clashes { policy, .. }
        | molframe::StructureRequest::BondLengthDeviations { policy, .. }
        | molframe::StructureRequest::CisPeptides { policy, .. }
        | molframe::StructureRequest::Planarity { policy, .. }
        | molframe::StructureRequest::Quality { policy }
        | molframe::StructureRequest::Valence { policy }
        | molframe::StructureRequest::Completeness { policy } => policy,
    }
}

pub(super) fn explain(py: Python<'_>, request: &molframe::StructureRequest) -> PyResult<Py<PyAny>> {
    let result = PyDict::new(py);
    result.set_item("operation", kind(request))?;
    result.set_item("backend", "rust")?;
    result.set_item("execution", "native")?;
    result.set_item("materializes", false)?;
    result.set_item(
        "policy_fingerprint",
        policy(request).fingerprint().to_string(),
    )?;
    result.set_item("complexity", complexity(request))?;
    Ok(result.unbind().into_any())
}

fn complexity(request: &molframe::StructureRequest) -> &'static str {
    match request {
        molframe::StructureRequest::PiStacking { .. }
        | molframe::StructureRequest::CationPi { .. }
        | molframe::StructureRequest::SecondaryStructure { .. } => "O(sparse²)",
        molframe::StructureRequest::GaussianNetworkModel { .. } => {
            "O(sites + local_neighbours + sites³)"
        }
        molframe::StructureRequest::WaterBridges { .. } => "O(hydrogen_bonds + bridges²)",
        molframe::StructureRequest::BasePairs { .. }
        | molframe::StructureRequest::HydrogenBonds { .. }
        | molframe::StructureRequest::SaltBridges { .. }
        | molframe::StructureRequest::ContactMap { .. }
        | molframe::StructureRequest::ChainInterface { .. }
        | molframe::StructureRequest::HalfSphereExposure { .. }
        | molframe::StructureRequest::NucleicTorsions { .. }
        | molframe::StructureRequest::Clashes { .. } => "O(n + local_neighbours)",
        molframe::StructureRequest::BondLengthDeviations { .. }
        | molframe::StructureRequest::CisPeptides { .. }
        | molframe::StructureRequest::Planarity { .. }
        | molframe::StructureRequest::Quality { .. }
        | molframe::StructureRequest::Valence { .. }
        | molframe::StructureRequest::Completeness { .. } => "O(n)",
    }
}

pub(super) fn to_dict(py: Python<'_>, request: &molframe::StructureRequest) -> PyResult<Py<PyAny>> {
    let result = PyDict::new(py);
    result.set_item("operation", kind(request))?;
    result.set_item(
        "policy",
        Py::new(py, PyAnalysisPolicy::from(policy(request).clone()))?,
    )?;
    super::config_fields::set_fields(py, &result, request)?;
    Ok(result.unbind().into_any())
}

pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<molframe::StructureRequest> {
    from_dict_with_provider(config, None)
}

pub(super) fn from_dict_with_provider(
    config: &Bound<'_, PyDict>,
    provider: Option<Arc<dyn molframe::ComponentProvider>>,
) -> PyResult<molframe::StructureRequest> {
    let operation: String = required(config, "operation")?.extract()?;
    let policy = match optional(config, "policy")? {
        Some(value) => value.extract::<PyAnalysisPolicy>()?.inner,
        None => molframe::AnalysisPolicy::default(),
    };
    if let Some(request) = from_analysis_dict(config, operation.as_str(), policy.clone(), provider)?
    {
        return Ok(request);
    }
    if let Some(request) = from_validation_dict(config, operation.as_str(), policy)? {
        return Ok(request);
    }
    Err(PyValueError::new_err(format!(
        "unknown structure operation {operation:?}"
    )))
}

fn from_analysis_dict(
    config: &Bound<'_, PyDict>,
    operation: &str,
    policy: molframe::AnalysisPolicy,
    provider: Option<Arc<dyn molframe::ComponentProvider>>,
) -> PyResult<Option<molframe::StructureRequest>> {
    match operation {
        "base_pairs" => {
            let Some(provider) = provider else {
                return Err(PyValueError::new_err(
                    "base_pairs requires a ComponentDictionary in the operation configuration",
                ));
            };
            Ok(Some(molframe::StructureRequest::BasePairs {
                provider,
                options: required(config, "options")?
                    .extract::<PyBasePairOptions>()?
                    .native(),
                policy,
            }))
        }
        "hydrogen_bonds" => Ok(Some(molframe::StructureRequest::HydrogenBonds {
            options: molframe::analysis::HydrogenBondOptions {
                maximum_donor_acceptor_distance: optional(config, "maximum_distance")?
                    .map_or(Ok(3.5), |value| value.extract())?,
                minimum_angle_degrees: optional(config, "minimum_angle")?
                    .map_or(Ok(120.0), |value| value.extract())?,
                backend: backend(config)?,
                periodic: optional(config, "periodic")?
                    .map_or(Ok(false), |value| value.extract())?,
            },
            policy,
        })),
        "salt_bridges" => Ok(Some(molframe::StructureRequest::SaltBridges {
            maximum_distance: optional(config, "maximum_distance")?
                .map_or(Ok(4.0), |value| value.extract())?,
            backend: backend(config)?,
            policy,
        })),
        "pi_stacking" => Ok(Some(molframe::StructureRequest::PiStacking {
            options: required(config, "options")?
                .extract::<PyPiStackingOptions>()?
                .native(),
            policy,
        })),
        "cation_pi" => Ok(Some(molframe::StructureRequest::CationPi {
            options: required(config, "options")?
                .extract::<PyCationPiOptions>()?
                .native(),
            policy,
        })),
        "water_bridges" => Ok(Some(molframe::StructureRequest::WaterBridges {
            options: required(config, "options")?
                .extract::<PyWaterBridgeOptions>()?
                .native(),
            policy,
        })),
        "contact_map" => Ok(Some(molframe::StructureRequest::ContactMap {
            cutoff: optional(config, "cutoff")?.map_or(Ok(4.5), |value| value.extract())?,
            minimum_separation: optional(config, "minimum_separation")?
                .map_or(Ok(0), |value| value.extract())?,
            backend: backend(config)?,
            policy,
        })),
        "chain_interface" => Ok(Some(molframe::StructureRequest::ChainInterface {
            first_chain: required(config, "first_chain")?
                .extract::<String>()?
                .into_boxed_str(),
            second_chain: required(config, "second_chain")?
                .extract::<String>()?
                .into_boxed_str(),
            cutoff: optional(config, "cutoff")?.map_or(Ok(4.5), |value| value.extract())?,
            backend: backend(config)?,
            policy,
        })),
        "secondary_structure" => Ok(Some(molframe::StructureRequest::SecondaryStructure {
            options: required(config, "options")?
                .extract::<PyDsspOptions>()?
                .native(),
            policy,
        })),
        "half_sphere_exposure" => Ok(Some(molframe::StructureRequest::HalfSphereExposure {
            radius: required(config, "radius")?.extract()?,
            backend: backend(config)?,
            policy,
        })),
        "gaussian_network_model" => Ok(Some(molframe::StructureRequest::GaussianNetworkModel {
            sites: required(config, "sites")?.extract::<PySelection>()?.inner,
            options: required(config, "options")?
                .extract::<PyGnmOptions>()?
                .native(),
            periodic: optional(config, "periodic")?.map_or(Ok(false), |value| value.extract())?,
            policy,
        })),
        "nucleic_torsions" => Ok(Some(molframe::StructureRequest::NucleicTorsions { policy })),
        "quality" => Ok(Some(molframe::StructureRequest::Quality { policy })),
        "valence" => Ok(Some(molframe::StructureRequest::Valence { policy })),
        "completeness" => Ok(Some(molframe::StructureRequest::Completeness { policy })),
        _ => Ok(None),
    }
}

fn from_validation_dict(
    config: &Bound<'_, PyDict>,
    operation: &str,
    policy: molframe::AnalysisPolicy,
) -> PyResult<Option<molframe::StructureRequest>> {
    match operation {
        "clashes" => Ok(Some(molframe::StructureRequest::Clashes {
            tolerance: optional(config, "tolerance")?.map_or(Ok(0.4), |value| value.extract())?,
            radii: radius_set(config)?,
            backend: backend(config)?,
            policy,
        })),
        "bond_length_deviations" => Ok(Some(molframe::StructureRequest::BondLengthDeviations {
            tolerance: required(config, "tolerance")?.extract()?,
            policy,
        })),
        "cis_peptides" => Ok(Some(molframe::StructureRequest::CisPeptides {
            threshold_degrees: required(config, "threshold_degrees")?.extract()?,
            policy,
        })),
        "planarity" => Ok(Some(molframe::StructureRequest::Planarity {
            options: required(config, "options")?
                .extract::<PyPlanarityOptions>()?
                .0,
            policy,
        })),
        _ => Ok(None),
    }
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

fn backend(config: &Bound<'_, PyDict>) -> PyResult<molframe::SpatialBackend> {
    let Some(value) = optional(config, "backend")? else {
        return Ok(molframe::SpatialBackend::Auto);
    };
    let value: PySpatialBackend = value.extract()?;
    Ok(value.into())
}

fn radius_set(config: &Bound<'_, PyDict>) -> PyResult<molframe::RadiusSet> {
    let Some(value) = optional(config, "radii")? else {
        return Ok(molframe::RadiusSet::Bondi);
    };
    let value: PyRadiusSet = value.extract()?;
    Ok(value.into())
}
