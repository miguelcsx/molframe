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

pub(super) fn kind(request: &pdbiox::StructureRequest) -> &'static str {
    match request {
        pdbiox::StructureRequest::BasePairs { .. } => "base_pairs",
        pdbiox::StructureRequest::HydrogenBonds { .. } => "hydrogen_bonds",
        pdbiox::StructureRequest::SaltBridges { .. } => "salt_bridges",
        pdbiox::StructureRequest::PiStacking { .. } => "pi_stacking",
        pdbiox::StructureRequest::CationPi { .. } => "cation_pi",
        pdbiox::StructureRequest::WaterBridges { .. } => "water_bridges",
        pdbiox::StructureRequest::ContactMap { .. } => "contact_map",
        pdbiox::StructureRequest::ChainInterface { .. } => "chain_interface",
        pdbiox::StructureRequest::SecondaryStructure { .. } => "secondary_structure",
        pdbiox::StructureRequest::HalfSphereExposure { .. } => "half_sphere_exposure",
        pdbiox::StructureRequest::NucleicTorsions { .. } => "nucleic_torsions",
        pdbiox::StructureRequest::GaussianNetworkModel { .. } => "gaussian_network_model",
        pdbiox::StructureRequest::Clashes { .. } => "clashes",
        pdbiox::StructureRequest::BondLengthDeviations { .. } => "bond_length_deviations",
        pdbiox::StructureRequest::CisPeptides { .. } => "cis_peptides",
        pdbiox::StructureRequest::Planarity { .. } => "planarity",
        pdbiox::StructureRequest::Quality { .. } => "quality",
        pdbiox::StructureRequest::Valence { .. } => "valence",
        pdbiox::StructureRequest::Completeness { .. } => "completeness",
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

fn policy(request: &pdbiox::StructureRequest) -> &pdbiox::AnalysisPolicy {
    match request {
        pdbiox::StructureRequest::BasePairs { policy, .. }
        | pdbiox::StructureRequest::HydrogenBonds { policy, .. }
        | pdbiox::StructureRequest::SaltBridges { policy, .. }
        | pdbiox::StructureRequest::PiStacking { policy, .. }
        | pdbiox::StructureRequest::CationPi { policy, .. }
        | pdbiox::StructureRequest::WaterBridges { policy, .. }
        | pdbiox::StructureRequest::ContactMap { policy, .. }
        | pdbiox::StructureRequest::ChainInterface { policy, .. }
        | pdbiox::StructureRequest::SecondaryStructure { policy, .. }
        | pdbiox::StructureRequest::HalfSphereExposure { policy, .. }
        | pdbiox::StructureRequest::NucleicTorsions { policy }
        | pdbiox::StructureRequest::GaussianNetworkModel { policy, .. }
        | pdbiox::StructureRequest::Clashes { policy, .. }
        | pdbiox::StructureRequest::BondLengthDeviations { policy, .. }
        | pdbiox::StructureRequest::CisPeptides { policy, .. }
        | pdbiox::StructureRequest::Planarity { policy, .. }
        | pdbiox::StructureRequest::Quality { policy }
        | pdbiox::StructureRequest::Valence { policy }
        | pdbiox::StructureRequest::Completeness { policy } => policy,
    }
}

pub(super) fn explain(py: Python<'_>, request: &pdbiox::StructureRequest) -> PyResult<Py<PyAny>> {
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

fn complexity(request: &pdbiox::StructureRequest) -> &'static str {
    match request {
        pdbiox::StructureRequest::PiStacking { .. }
        | pdbiox::StructureRequest::CationPi { .. }
        | pdbiox::StructureRequest::SecondaryStructure { .. } => "O(sparse²)",
        pdbiox::StructureRequest::GaussianNetworkModel { .. } => {
            "O(sites + local_neighbours + sites³)"
        }
        pdbiox::StructureRequest::WaterBridges { .. } => "O(hydrogen_bonds + bridges²)",
        pdbiox::StructureRequest::BasePairs { .. }
        | pdbiox::StructureRequest::HydrogenBonds { .. }
        | pdbiox::StructureRequest::SaltBridges { .. }
        | pdbiox::StructureRequest::ContactMap { .. }
        | pdbiox::StructureRequest::ChainInterface { .. }
        | pdbiox::StructureRequest::HalfSphereExposure { .. }
        | pdbiox::StructureRequest::NucleicTorsions { .. }
        | pdbiox::StructureRequest::Clashes { .. } => "O(n + local_neighbours)",
        pdbiox::StructureRequest::BondLengthDeviations { .. }
        | pdbiox::StructureRequest::CisPeptides { .. }
        | pdbiox::StructureRequest::Planarity { .. }
        | pdbiox::StructureRequest::Quality { .. }
        | pdbiox::StructureRequest::Valence { .. }
        | pdbiox::StructureRequest::Completeness { .. } => "O(n)",
    }
}

pub(super) fn to_dict(py: Python<'_>, request: &pdbiox::StructureRequest) -> PyResult<Py<PyAny>> {
    let result = PyDict::new(py);
    result.set_item("operation", kind(request))?;
    result.set_item(
        "policy",
        Py::new(py, PyAnalysisPolicy::from(policy(request).clone()))?,
    )?;
    super::config_fields::set_fields(py, &result, request)?;
    Ok(result.unbind().into_any())
}

pub(super) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<pdbiox::StructureRequest> {
    from_dict_with_provider(config, None)
}

pub(super) fn from_dict_with_provider(
    config: &Bound<'_, PyDict>,
    provider: Option<Arc<dyn pdbiox::ComponentProvider>>,
) -> PyResult<pdbiox::StructureRequest> {
    let operation: String = required(config, "operation")?.extract()?;
    let policy = match optional(config, "policy")? {
        Some(value) => value.extract::<PyAnalysisPolicy>()?.inner,
        None => pdbiox::AnalysisPolicy::default(),
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
    policy: pdbiox::AnalysisPolicy,
    provider: Option<Arc<dyn pdbiox::ComponentProvider>>,
) -> PyResult<Option<pdbiox::StructureRequest>> {
    match operation {
        "base_pairs" => {
            let Some(provider) = provider else {
                return Err(PyValueError::new_err(
                    "base_pairs requires a ComponentDictionary in the operation configuration",
                ));
            };
            Ok(Some(pdbiox::StructureRequest::BasePairs {
                provider,
                options: required(config, "options")?
                    .extract::<PyBasePairOptions>()?
                    .native(),
                policy,
            }))
        }
        "hydrogen_bonds" => Ok(Some(pdbiox::StructureRequest::HydrogenBonds {
            options: pdbiox::analysis::HydrogenBondOptions {
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
        "salt_bridges" => Ok(Some(pdbiox::StructureRequest::SaltBridges {
            maximum_distance: optional(config, "maximum_distance")?
                .map_or(Ok(4.0), |value| value.extract())?,
            backend: backend(config)?,
            policy,
        })),
        "pi_stacking" => Ok(Some(pdbiox::StructureRequest::PiStacking {
            options: required(config, "options")?
                .extract::<PyPiStackingOptions>()?
                .native(),
            policy,
        })),
        "cation_pi" => Ok(Some(pdbiox::StructureRequest::CationPi {
            options: required(config, "options")?
                .extract::<PyCationPiOptions>()?
                .native(),
            policy,
        })),
        "water_bridges" => Ok(Some(pdbiox::StructureRequest::WaterBridges {
            options: required(config, "options")?
                .extract::<PyWaterBridgeOptions>()?
                .native(),
            policy,
        })),
        "contact_map" => Ok(Some(pdbiox::StructureRequest::ContactMap {
            cutoff: optional(config, "cutoff")?.map_or(Ok(4.5), |value| value.extract())?,
            minimum_separation: optional(config, "minimum_separation")?
                .map_or(Ok(0), |value| value.extract())?,
            backend: backend(config)?,
            policy,
        })),
        "chain_interface" => Ok(Some(pdbiox::StructureRequest::ChainInterface {
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
        "secondary_structure" => Ok(Some(pdbiox::StructureRequest::SecondaryStructure {
            options: required(config, "options")?
                .extract::<PyDsspOptions>()?
                .native(),
            policy,
        })),
        "half_sphere_exposure" => Ok(Some(pdbiox::StructureRequest::HalfSphereExposure {
            radius: required(config, "radius")?.extract()?,
            backend: backend(config)?,
            policy,
        })),
        "gaussian_network_model" => Ok(Some(pdbiox::StructureRequest::GaussianNetworkModel {
            sites: required(config, "sites")?.extract::<PySelection>()?.inner,
            options: required(config, "options")?
                .extract::<PyGnmOptions>()?
                .native(),
            periodic: optional(config, "periodic")?.map_or(Ok(false), |value| value.extract())?,
            policy,
        })),
        "nucleic_torsions" => Ok(Some(pdbiox::StructureRequest::NucleicTorsions { policy })),
        "quality" => Ok(Some(pdbiox::StructureRequest::Quality { policy })),
        "valence" => Ok(Some(pdbiox::StructureRequest::Valence { policy })),
        "completeness" => Ok(Some(pdbiox::StructureRequest::Completeness { policy })),
        _ => Ok(None),
    }
}

fn from_validation_dict(
    config: &Bound<'_, PyDict>,
    operation: &str,
    policy: pdbiox::AnalysisPolicy,
) -> PyResult<Option<pdbiox::StructureRequest>> {
    match operation {
        "clashes" => Ok(Some(pdbiox::StructureRequest::Clashes {
            tolerance: optional(config, "tolerance")?.map_or(Ok(0.4), |value| value.extract())?,
            radii: radius_set(config)?,
            backend: backend(config)?,
            policy,
        })),
        "bond_length_deviations" => Ok(Some(pdbiox::StructureRequest::BondLengthDeviations {
            tolerance: required(config, "tolerance")?.extract()?,
            policy,
        })),
        "cis_peptides" => Ok(Some(pdbiox::StructureRequest::CisPeptides {
            threshold_degrees: required(config, "threshold_degrees")?.extract()?,
            policy,
        })),
        "planarity" => Ok(Some(pdbiox::StructureRequest::Planarity {
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

fn backend(config: &Bound<'_, PyDict>) -> PyResult<pdbiox::SpatialBackend> {
    let Some(value) = optional(config, "backend")? else {
        return Ok(pdbiox::SpatialBackend::Auto);
    };
    let value: PySpatialBackend = value.extract()?;
    Ok(value.into())
}

fn radius_set(config: &Bound<'_, PyDict>) -> PyResult<pdbiox::RadiusSet> {
    let Some(value) = optional(config, "radii")? else {
        return Ok(pdbiox::RadiusSet::Bondi);
    };
    let value: PyRadiusSet = value.extract()?;
    Ok(value.into())
}
