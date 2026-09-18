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
    request: &molframe::StructureRequest,
) -> PyResult<()> {
    match request {
        molframe::StructureRequest::BasePairs { .. }
        | molframe::StructureRequest::HydrogenBonds { .. }
        | molframe::StructureRequest::SaltBridges { .. }
        | molframe::StructureRequest::PiStacking { .. }
        | molframe::StructureRequest::CationPi { .. }
        | molframe::StructureRequest::WaterBridges { .. }
        | molframe::StructureRequest::ContactMap { .. }
        | molframe::StructureRequest::ChainInterface { .. } => {
            set_interaction_fields(py, result, request)
        }
        molframe::StructureRequest::SecondaryStructure { .. }
        | molframe::StructureRequest::HalfSphereExposure { .. }
        | molframe::StructureRequest::GaussianNetworkModel { .. }
        | molframe::StructureRequest::NucleicTorsions { .. }
        | molframe::StructureRequest::Quality { .. }
        | molframe::StructureRequest::Valence { .. }
        | molframe::StructureRequest::Completeness { .. } => {
            set_structure_fields(py, result, request)
        }
        molframe::StructureRequest::Clashes { .. }
        | molframe::StructureRequest::BondLengthDeviations { .. }
        | molframe::StructureRequest::CisPeptides { .. }
        | molframe::StructureRequest::Planarity { .. } => {
            set_validation_fields(py, result, request)
        }
    }
}

fn set_interaction_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    request: &molframe::StructureRequest,
) -> PyResult<()> {
    match request {
        molframe::StructureRequest::BasePairs { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyBasePairOptions::from_native(*options))?,
            )?;
        }
        molframe::StructureRequest::HydrogenBonds { options, .. } => {
            result.set_item("maximum_distance", options.maximum_donor_acceptor_distance)?;
            result.set_item("minimum_angle", options.minimum_angle_degrees)?;
            result.set_item("backend", PySpatialBackend::from(options.backend))?;
            result.set_item("periodic", options.periodic)?;
        }
        molframe::StructureRequest::SaltBridges {
            maximum_distance,
            backend,
            ..
        } => {
            result.set_item("maximum_distance", maximum_distance)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        molframe::StructureRequest::PiStacking { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyPiStackingOptions::from_native_parts(*options))?,
            )?;
        }
        molframe::StructureRequest::CationPi { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyCationPiOptions::from_native_parts(*options))?,
            )?;
        }
        molframe::StructureRequest::WaterBridges { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyWaterBridgeOptions::from_native_parts(*options))?,
            )?;
        }
        molframe::StructureRequest::ContactMap {
            cutoff,
            minimum_separation,
            backend,
            ..
        } => {
            result.set_item("cutoff", cutoff)?;
            result.set_item("minimum_separation", minimum_separation)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        molframe::StructureRequest::ChainInterface {
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
    request: &molframe::StructureRequest,
) -> PyResult<()> {
    match request {
        molframe::StructureRequest::SecondaryStructure { options, .. } => {
            result.set_item(
                "options",
                Py::new(py, PyDsspOptions::from_native_parts(options.clone()))?,
            )?;
        }
        molframe::StructureRequest::HalfSphereExposure {
            radius, backend, ..
        } => {
            result.set_item("radius", radius)?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        molframe::StructureRequest::GaussianNetworkModel {
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
        molframe::StructureRequest::NucleicTorsions { .. }
        | molframe::StructureRequest::Quality { .. }
        | molframe::StructureRequest::Valence { .. }
        | molframe::StructureRequest::Completeness { .. } => {}
        _ => Err(PyValueError::new_err(
            "structure analysis field projection mismatch",
        ))?,
    }
    Ok(())
}

fn set_validation_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    request: &molframe::StructureRequest,
) -> PyResult<()> {
    match request {
        molframe::StructureRequest::Clashes {
            tolerance,
            radii,
            backend,
            ..
        } => {
            result.set_item("tolerance", tolerance)?;
            result.set_item("radii", PyRadiusSet::from(*radii))?;
            result.set_item("backend", PySpatialBackend::from(*backend))?;
        }
        molframe::StructureRequest::BondLengthDeviations { tolerance, .. } => {
            result.set_item("tolerance", tolerance)?;
        }
        molframe::StructureRequest::CisPeptides {
            threshold_degrees, ..
        } => {
            result.set_item("threshold_degrees", threshold_degrees)?;
        }
        molframe::StructureRequest::Planarity { options, .. } => {
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
