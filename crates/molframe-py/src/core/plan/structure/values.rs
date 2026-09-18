//! Typed projections for structure-plan results.

use crate::analysis::PyBondDeviation;
use crate::analysis::{
    PyBasePair, PyCationPi, PyChainCompleteness, PyCisPeptide, PyClash, PyContactMap,
    PyGaussianNetworkModel, PyHalfSphereExposure, PyHydrogenBond, PyNucleicTorsions, PyPiStacking,
    PyPlanarityFlag, PyQualityFlag, PySaltBridge, PySseRecord, PyValenceError, PyWaterBridge,
};
use crate::contract::{PyAnalysis, analysis_with_value};
use pyo3::prelude::*;
use pyo3::types::PyList;

macro_rules! list_analysis {
    ($function:ident, $rust:ty, $python:ty) => {
        fn $function(
            py: Python<'_>,
            analysis: molframe::Analysis<Vec<$rust>>,
        ) -> PyResult<PyAnalysis> {
            analysis_with_value(py, analysis, |py, values| {
                let items = values
                    .into_iter()
                    .map(|value| Py::new(py, <$python>::from(value)))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok(PyList::new(py, items)?.unbind().into_any())
            })
        }
    };
}

list_analysis!(
    hydrogen_bonds,
    molframe::analysis::HydrogenBond,
    PyHydrogenBond
);
list_analysis!(base_pairs, molframe::analysis::BasePair, PyBasePair);
list_analysis!(salt_bridges, molframe::analysis::SaltBridge, PySaltBridge);
list_analysis!(pi_stacking, molframe::analysis::PiStacking, PyPiStacking);
list_analysis!(cation_pi, molframe::analysis::CationPi, PyCationPi);
list_analysis!(
    water_bridges,
    molframe::analysis::WaterBridge,
    PyWaterBridge
);
list_analysis!(
    secondary_structure,
    molframe::analysis::SseRecord,
    PySseRecord
);
list_analysis!(
    half_sphere_exposure,
    molframe::analysis::HalfSphereExposure,
    PyHalfSphereExposure
);
list_analysis!(
    nucleic_torsions,
    molframe::analysis::NucleicTorsions,
    PyNucleicTorsions
);

fn gaussian_network_model(
    py: Python<'_>,
    analysis: molframe::Analysis<molframe::analysis::GaussianNetworkModel>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, value| {
        let model = PyGaussianNetworkModel::new(py, value)?;
        Ok(Py::new(py, model)?.into_any())
    })
}

list_analysis!(clashes, molframe::validate::Clash, PyClash);
list_analysis!(
    bond_length_deviations,
    molframe::validate::BondDeviation,
    PyBondDeviation
);
list_analysis!(cis_peptides, molframe::validate::CisPeptide, PyCisPeptide);
list_analysis!(
    planarity,
    molframe::validate::PlanarityFlag,
    PyPlanarityFlag
);
list_analysis!(quality, molframe::validate::QualityFlag, PyQualityFlag);
list_analysis!(valence, molframe::validate::ValenceError, PyValenceError);

fn contact_map(
    py: Python<'_>,
    analysis: molframe::Analysis<molframe::analysis::ContactMap>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, value| {
        Ok(Py::new(py, PyContactMap::from(value))?.into_any())
    })
}

fn chain_interface(
    py: Python<'_>,
    analysis: molframe::Analysis<Vec<molframe::ResidueIndex>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        Ok(
            PyList::new(py, values.into_iter().map(molframe::ResidueIndex::get))?
                .unbind()
                .into_any(),
        )
    })
}

fn completeness(
    py: Python<'_>,
    structure: Option<&molframe::Structure>,
    analysis: molframe::Analysis<Vec<molframe::validate::ChainCompleteness>>,
) -> PyResult<PyAnalysis> {
    let Some(structure) = structure else {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "completeness results require the retained structure snapshot",
        ));
    };
    let projected = analysis
        .value
        .into_iter()
        .map(|value| crate::analysis::project_completeness(structure, value))
        .collect::<Result<Vec<PyChainCompleteness>, _>>()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    let analysis = molframe::Analysis {
        value: projected,
        status: analysis.status,
        coverage: analysis.coverage,
        warnings: analysis.warnings,
        assumptions: analysis.assumptions,
        provenance: analysis.provenance,
    };
    analysis_with_value(py, analysis, |py, values| {
        let items = values
            .into_iter()
            .map(|value| Py::new(py, value))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, items)?.unbind().into_any())
    })
}

pub(super) fn to_python(
    py: Python<'_>,
    value: molframe::StructureValue,
    structure: Option<&molframe::Structure>,
) -> PyResult<PyAnalysis> {
    match value {
        molframe::StructureValue::BasePairs(value) => base_pairs(py, value),
        molframe::StructureValue::HydrogenBonds(value) => hydrogen_bonds(py, value),
        molframe::StructureValue::SaltBridges(value) => salt_bridges(py, value),
        molframe::StructureValue::PiStacking(value) => pi_stacking(py, value),
        molframe::StructureValue::CationPi(value) => cation_pi(py, value),
        molframe::StructureValue::WaterBridges(value) => water_bridges(py, value),
        molframe::StructureValue::ContactMap(value) => contact_map(py, value),
        molframe::StructureValue::ChainInterface(value) => chain_interface(py, value),
        molframe::StructureValue::SecondaryStructure(value) => secondary_structure(py, value),
        molframe::StructureValue::HalfSphereExposure(value) => half_sphere_exposure(py, value),
        molframe::StructureValue::NucleicTorsions(value) => nucleic_torsions(py, value),
        molframe::StructureValue::GaussianNetworkModel(value) => gaussian_network_model(py, value),
        molframe::StructureValue::Clashes(value) => clashes(py, value),
        molframe::StructureValue::BondLengthDeviations(value) => bond_length_deviations(py, value),
        molframe::StructureValue::CisPeptides(value) => cis_peptides(py, value),
        molframe::StructureValue::Planarity(value) => planarity(py, value),
        molframe::StructureValue::Quality(value) => quality(py, value),
        molframe::StructureValue::Valence(value) => valence(py, value),
        molframe::StructureValue::Completeness(value) => completeness(py, structure, value),
    }
}
