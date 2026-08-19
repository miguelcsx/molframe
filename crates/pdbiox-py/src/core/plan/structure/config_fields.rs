//! Python projections for typed structure-operation fields.

use crate::analysis::{
    PyBasePairOptions, PyCationPiOptions, PyDsspOptions, PyGnmOptions, PyPiStackingOptions,
    PyPlanarityOptions, PyWaterBridgeOptions,
};
use crate::chemistry::PyRadiusSet;
use crate::graph::PySpatialBackend;
use crate::query::PySelection;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

pub(super) fn set_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    request: &pdbiox::StructureRequest,
) -> PyResult<()> {
    match request {
        pdbiox::StructureRequest::BasePairs { .. }
        | pdbiox::StructureRequest::HydrogenBonds { .. }
        | pdbiox::StructureRequest::SaltBridges { .. }
        | pdbiox::StructureRequest::PiStacking { .. }
        | pdbiox::StructureRequest::CationPi { .. }
        | pdbiox::StructureRequest::WaterBridges { .. }
        | pdbiox::StructureRequest::ContactMap { .. }
        | pdbiox::StructureRequest::ChainInterface { .. } => {
            set_interaction_fields(py, result, request)
        }
        pdbiox::StructureRequest::SecondaryStructure { .. }
        | pdbiox::StructureRequest::HalfSphereExposure { .. }
        | pdbiox::StructureRequest::GaussianNetworkModel { .. }
        | pdbiox::StructureRequest::NucleicTorsions { .. }
        | pdbiox::StructureRequest::Quality { .. }
        | pdbiox::StructureRequest::Valence { .. }
        | pdbiox::StructureRequest::Completeness { .. } => {
            set_structure_fields(py, result, request)
        }
        pdbiox::StructureRequest::Clashes { .. }
        | pdbiox::StructureRequest::BondLengthDeviations { .. }
        | pdbiox::StructureRequest::CisPeptides { .. }
        | pdbiox::StructureRequest::Planarity { .. } => set_validation_fields(py, result, request),
    }
}

fn set_interaction_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    request: &pdbiox::StructureRequest,
) -> PyResult<()> {
    match request {
        pdbiox::StructureRequest::BasePairs { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyBasePairOptions::from_native(*options))?,
            )?;
        }
        pdbiox::StructureRequest::HydrogenBonds { options, .. } => {
            result.set_item("maximum_distance", options.maximum_donor_acceptor_distance)?;
            result.set_item("minimum_angle", options.minimum_angle_degrees)?;
            result.set_item("backend", PySpatialBackend::from(options.backend))?;
            result.set_item("periodic", options.periodic)?;
        }
        pdbiox::StructureRequest::SaltBridges {
            maximum_distance,
            backend,
            ..
        } => {
            result.set_item("maximum_distance", maximum_distance)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        pdbiox::StructureRequest::PiStacking { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyPiStackingOptions::from_native_parts(*options))?,
            )?;
        }
        pdbiox::StructureRequest::CationPi { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyCationPiOptions::from_native_parts(*options))?,
            )?;
        }
        pdbiox::StructureRequest::WaterBridges { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyWaterBridgeOptions::from_native_parts(*options))?,
            )?;
        }
        pdbiox::StructureRequest::ContactMap {
            cutoff,
            minimum_separation,
            backend,
            ..
        } => {
            result.set_item("cutoff", cutoff)?;
            result.set_item("minimum_separation", minimum_separation)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        pdbiox::StructureRequest::ChainInterface {
            first_chain,
            second_chain,
            cutoff,
            backend,
            ..
        } => {
            result.set_item("first_chain", first_chain.as_ref())?;
            result.set_item("second_chain", second_chain.as_ref())?;
            result.set_item("cutoff", cutoff)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        _ => Err(PyValueError::new_err(
            "structure interaction field projection mismatch",
        ))?,
    }
    Ok(())
}

fn set_structure_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    request: &pdbiox::StructureRequest,
) -> PyResult<()> {
    match request {
        pdbiox::StructureRequest::SecondaryStructure { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyDsspOptions::from_native_parts(options.clone()))?,
            )?;
        }
        pdbiox::StructureRequest::HalfSphereExposure {
            radius, backend, ..
        } => {
            result.set_item("radius", radius)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        pdbiox::StructureRequest::GaussianNetworkModel {
            sites,
            options,
            periodic,
            ..
        } => {
            result.set_item(
                "sites",
                Py::new(
                    py,
                    PySelection {
                        inner: sites.clone(),
                    },
                )?,
            )?;
            result.set_item("options", Py::new(py, PyGnmOptions::from_native(*options))?)?;
            result.set_item("periodic", periodic)?;
        }
        pdbiox::StructureRequest::NucleicTorsions { .. }
        | pdbiox::StructureRequest::Quality { .. }
        | pdbiox::StructureRequest::Valence { .. }
        | pdbiox::StructureRequest::Completeness { .. } => {}
        _ => Err(PyValueError::new_err(
            "structure analysis field projection mismatch",
        ))?,
    }
    Ok(())
}

fn set_validation_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    request: &pdbiox::StructureRequest,
) -> PyResult<()> {
    match request {
        pdbiox::StructureRequest::Clashes {
            tolerance,
            radii,
            backend,
            ..
        } => {
            result.set_item("tolerance", tolerance)?;
            result.set_item("radii", PyRadiusSet::from(*radii))?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        pdbiox::StructureRequest::BondLengthDeviations { tolerance, .. } => {
            result.set_item("tolerance", tolerance)?;
        }
        pdbiox::StructureRequest::CisPeptides {
            threshold_degrees, ..
        } => {
            result.set_item("threshold_degrees", threshold_degrees)?;
        }
        pdbiox::StructureRequest::Planarity { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyPlanarityOptions::from_native_parts(*options))?,
            )?;
        }
        _ => Err(PyValueError::new_err(
            "structure validation field projection mismatch",
        ))?,
    }
    Ok(())
}
