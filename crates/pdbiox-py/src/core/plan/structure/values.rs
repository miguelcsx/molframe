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
            analysis: pdbiox::Analysis<Vec<$rust>>,
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
    pdbiox::analysis::HydrogenBond,
    PyHydrogenBond
);
list_analysis!(base_pairs, pdbiox::analysis::BasePair, PyBasePair);
list_analysis!(salt_bridges, pdbiox::analysis::SaltBridge, PySaltBridge);
list_analysis!(pi_stacking, pdbiox::analysis::PiStacking, PyPiStacking);
list_analysis!(cation_pi, pdbiox::analysis::CationPi, PyCationPi);
list_analysis!(water_bridges, pdbiox::analysis::WaterBridge, PyWaterBridge);
list_analysis!(
    secondary_structure,
    pdbiox::analysis::SseRecord,
    PySseRecord
);
list_analysis!(
    half_sphere_exposure,
    pdbiox::analysis::HalfSphereExposure,
    PyHalfSphereExposure
);
list_analysis!(
    nucleic_torsions,
    pdbiox::analysis::NucleicTorsions,
    PyNucleicTorsions
);

fn gaussian_network_model(
    py: Python<'_>,
    analysis: pdbiox::Analysis<pdbiox::analysis::GaussianNetworkModel>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, value| {
        let model = PyGaussianNetworkModel::new(py, value)?;
        Ok(Py::new(py, model)?.into_any())
    })
}

list_analysis!(clashes, pdbiox::validate::Clash, PyClash);
list_analysis!(
    bond_length_deviations,
    pdbiox::validate::BondDeviation,
    PyBondDeviation
);
list_analysis!(cis_peptides, pdbiox::validate::CisPeptide, PyCisPeptide);
list_analysis!(planarity, pdbiox::validate::PlanarityFlag, PyPlanarityFlag);
list_analysis!(quality, pdbiox::validate::QualityFlag, PyQualityFlag);
list_analysis!(valence, pdbiox::validate::ValenceError, PyValenceError);

fn contact_map(
    py: Python<'_>,
    analysis: pdbiox::Analysis<pdbiox::analysis::ContactMap>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, value| {
        Ok(Py::new(py, PyContactMap::from(value))?.into_any())
    })
}

fn chain_interface(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::ResidueIndex>>,
) -> PyResult<PyAnalysis> {
    analysis_with_value(py, analysis, |py, values| {
        Ok(
            PyList::new(py, values.into_iter().map(pdbiox::ResidueIndex::get))?
                .unbind()
                .into_any(),
        )
    })
}

fn completeness(
    py: Python<'_>,
    structure: Option<&pdbiox::Structure>,
    analysis: pdbiox::Analysis<Vec<pdbiox::validate::ChainCompleteness>>,
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
    let analysis = pdbiox::Analysis {
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
    value: pdbiox::StructureValue,
    structure: Option<&pdbiox::Structure>,
) -> PyResult<PyAnalysis> {
    match value {
        pdbiox::StructureValue::BasePairs(value) => base_pairs(py, value),
        pdbiox::StructureValue::HydrogenBonds(value) => hydrogen_bonds(py, value),
        pdbiox::StructureValue::SaltBridges(value) => salt_bridges(py, value),
        pdbiox::StructureValue::PiStacking(value) => pi_stacking(py, value),
        pdbiox::StructureValue::CationPi(value) => cation_pi(py, value),
        pdbiox::StructureValue::WaterBridges(value) => water_bridges(py, value),
        pdbiox::StructureValue::ContactMap(value) => contact_map(py, value),
        pdbiox::StructureValue::ChainInterface(value) => chain_interface(py, value),
        pdbiox::StructureValue::SecondaryStructure(value) => secondary_structure(py, value),
        pdbiox::StructureValue::HalfSphereExposure(value) => half_sphere_exposure(py, value),
        pdbiox::StructureValue::NucleicTorsions(value) => nucleic_torsions(py, value),
        pdbiox::StructureValue::GaussianNetworkModel(value) => gaussian_network_model(py, value),
        pdbiox::StructureValue::Clashes(value) => clashes(py, value),
        pdbiox::StructureValue::BondLengthDeviations(value) => bond_length_deviations(py, value),
        pdbiox::StructureValue::CisPeptides(value) => cis_peptides(py, value),
        pdbiox::StructureValue::Planarity(value) => planarity(py, value),
        pdbiox::StructureValue::Quality(value) => quality(py, value),
        pdbiox::StructureValue::Valence(value) => valence(py, value),
        pdbiox::StructureValue::Completeness(value) => completeness(py, structure, value),
    }
}
