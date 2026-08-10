//! Policy-bound structure analysis and validation entry points.

use super::{
    PyBondDeviation, PyChainCompleteness, PyCisPeptide, PyClash, PyContact, PyHalfSphereExposure,
    PyNucleicTorsions, PyPlanarityFlag, PyPlanarityOptions, PyQualityFlag, PyValenceError,
};
use crate::chemistry::PyRadiusSet;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::graph::PySpatialBackend;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

macro_rules! class_list {
    ($name:ident, $rust:ty, $python:ty) => {
        fn $name(py: Python<'_>, values: Vec<$rust>) -> PyResult<Py<PyAny>> {
            let items = values
                .into_iter()
                .map(|value| Py::new(py, <$python>::from(value)))
                .collect::<PyResult<Vec<_>>>()?;
            Ok(PyList::new(py, items)?.unbind().into_any())
        }
    };
}

class_list!(contact_list, pdbiox::analysis::Contact, PyContact);
class_list!(
    hse_list,
    pdbiox::analysis::HalfSphereExposure,
    PyHalfSphereExposure
);
class_list!(
    torsion_list,
    pdbiox::analysis::NucleicTorsions,
    PyNucleicTorsions
);
class_list!(clash_list, pdbiox::validate::Clash, PyClash);
class_list!(bond_list, pdbiox::validate::BondDeviation, PyBondDeviation);
class_list!(cis_list, pdbiox::validate::CisPeptide, PyCisPeptide);
class_list!(
    planarity_list,
    pdbiox::validate::PlanarityFlag,
    PyPlanarityFlag
);
class_list!(quality_list, pdbiox::validate::QualityFlag, PyQualityFlag);
class_list!(valence_list, pdbiox::validate::ValenceError, PyValenceError);
class_list!(completeness_list, PyChainCompleteness, PyChainCompleteness);

fn index_list(py: Python<'_>, values: Vec<pdbiox::ResidueIndex>) -> PyResult<Py<PyAny>> {
    Ok(
        PyList::new(py, values.into_iter().map(pdbiox::ResidueIndex::get))?
            .unbind()
            .into_any(),
    )
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

macro_rules! governed_list {
    ($name:ident, $converter:expr, $kernel:expr) => {
        #[pyfunction]
        pub(crate) fn $name(
            py: Python<'_>,
            structure: &PyStructure,
            policy: &PyAnalysisPolicy,
        ) -> PyResult<PyAnalysis> {
            let structure = structure.structure().clone();
            let policy = policy.inner.clone();
            let analysis = py
                .detach(|| {
                    let kernel = $kernel;
                    pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
                })
                .map_err(value_error)?;
            analysis_with_value(py, analysis, $converter)
        }
    };
}

#[pyfunction]
pub(crate) fn analyse_contacts(
    py: Python<'_>,
    structure: &PyStructure,
    cutoff: f32,
    backend: PySpatialBackend,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::analysis::contacts_kernel(cutoff, backend.into());
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, contact_list)
}

#[pyfunction]
pub(crate) fn analyse_chain_interface(
    py: Python<'_>,
    structure: &PyStructure,
    first_chain: &str,
    second_chain: &str,
    cutoff: f32,
    backend: PySpatialBackend,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::analysis::chain_interface_kernel(
                first_chain,
                second_chain,
                cutoff,
                backend.into(),
            );
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, index_list)
}

#[pyfunction]
pub(crate) fn analyse_half_sphere_exposure(
    py: Python<'_>,
    structure: &PyStructure,
    radius: f32,
    backend: PySpatialBackend,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::analysis::half_sphere_exposure_kernel(radius, backend.into());
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, hse_list)
}

governed_list!(
    analyse_nucleic_torsions,
    torsion_list,
    pdbiox::analysis::nucleic_torsions_kernel()
);

#[pyfunction]
pub(crate) fn validate_clashes(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    radii: PyRadiusSet,
    backend: PySpatialBackend,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::validate::clashes_kernel(tolerance, radii.into(), backend.into());
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, clash_list)
}

#[pyfunction]
pub(crate) fn validate_bond_lengths(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::validate::bond_length_deviations_kernel(tolerance);
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, bond_list)
}

#[pyfunction]
pub(crate) fn validate_cis_peptides(
    py: Python<'_>,
    structure: &PyStructure,
    threshold_degrees: f64,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::validate::cis_peptides_kernel(threshold_degrees);
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, cis_list)
}

#[pyfunction]
pub(crate) fn validate_planarity(
    py: Python<'_>,
    structure: &PyStructure,
    options: PyPlanarityOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::validate::planarity_kernel(options.0);
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, planarity_list)
}

governed_list!(
    validate_quality,
    quality_list,
    pdbiox::validate::quality_flags_kernel()
);
governed_list!(
    validate_valence,
    valence_list,
    pdbiox::validate::valence_kernel()
);

#[pyfunction]
pub(crate) fn validate_completeness(
    py: Python<'_>,
    structure: &PyStructure,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let retained = structure.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            let kernel = pdbiox::validate::completeness_kernel();
            pdbiox::analysis::analyse_structure(&structure, &policy, &kernel)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, values| {
        let projected = values
            .into_iter()
            .map(|value| super::general_validation::project_completeness(&retained, value))
            .collect::<Result<Vec<PyChainCompleteness>, _>>()
            .map_err(value_error)?;
        completeness_list(py, projected)
    })
}
