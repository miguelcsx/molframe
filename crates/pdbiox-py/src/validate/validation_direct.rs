//! Direct validation functions mirroring the facade's reusable kernels.

use super::{
    PyAltlocOccupancyOptions, PyBondDeviation, PyChainCompleteness, PyChiralityOptions,
    PyChiralityReport, PyCisPeptide, PyClash, PyNucleicGeometryPolicy, PyNucleicGeometryRecord,
    PyPlanarityFlag, PyPlanarityOptions, PyQualityFlag, PyRamachandranOptions,
    PyRamachandranRecord, PyReferenceGeometryOptions, PyReferenceGeometryReport,
    PyReferenceLibrary, PyRotamerOptions, PyRotamerProfile, PyRotamerReport, PyValenceError,
};
use crate::analysis::general_validation::project_completeness;
use crate::chemistry::PyComponentDictionary as RootComponentDictionary;
use crate::chemistry::{PyPolymerRoleProfile, PyRadiusSet};
use crate::graph::PySpatialBackend;
use crate::query::{PyAnalysisPolicy, PyNamespace};
use crate::structure::PyStructure;
use pyo3::prelude::*;

fn value_error(error: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}

#[pyfunction]
#[pyo3(signature = (structure, tolerance))]
pub(crate) fn bond_length_deviations(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
) -> Vec<PyBondDeviation> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::bond_length_deviations(&structure, tolerance))
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
#[pyo3(signature = (structure, tolerance))]
pub(crate) fn ligand_geometry_outliers(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
) -> Vec<PyBondDeviation> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::ligand_geometry_outliers(&structure, tolerance))
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
#[pyo3(signature = (structure, tolerance))]
pub(crate) fn ligand_geometry(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
) -> super::PyLigandGeometryReport {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::ligand_geometry(&structure, tolerance))
        .into()
}

#[pyfunction]
#[pyo3(signature = (structure, tolerance, radius_set, *, backend=PySpatialBackend::Auto))]
pub(crate) fn clashes(
    py: Python<'_>,
    structure: &PyStructure,
    tolerance: f32,
    radius_set: PyRadiusSet,
    backend: PySpatialBackend,
) -> PyResult<Vec<PyClash>> {
    let structure = structure.structure().clone();
    py.detach(move || {
        pdbiox::validate::clashes(&structure, tolerance, radius_set.into(), backend.into())
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn completeness(
    py: Python<'_>,
    structure: &PyStructure,
    namespace: PyNamespace,
) -> PyResult<Vec<PyChainCompleteness>> {
    let structure = structure.structure().clone();
    let retained = structure.clone();
    py.detach(move || pdbiox::validate::completeness(&structure, namespace.into()))
        .map_err(value_error)?
        .into_iter()
        .map(|value| project_completeness(&retained, value).map_err(value_error))
        .collect()
}

#[pyfunction]
pub(crate) fn cis_peptides(
    py: Python<'_>,
    structure: &PyStructure,
    threshold_degrees: f64,
) -> PyResult<Vec<PyCisPeptide>> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::cis_peptides(&structure, threshold_degrees))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn nonplanar_aromatic_rings(
    py: Python<'_>,
    structure: &PyStructure,
    options: PyPlanarityOptions,
) -> PyResult<Vec<PyPlanarityFlag>> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::nonplanar_aromatic_rings(&structure, options.0))
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn overvalent_atoms(py: Python<'_>, structure: &PyStructure) -> Vec<PyValenceError> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::overvalent_atoms(&structure))
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
pub(crate) fn quality_flags(py: Python<'_>, structure: &PyStructure) -> Vec<PyQualityFlag> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::validate::quality_flags(&structure))
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
#[pyo3(signature = (structure, options, *, outliers_only=false))]
pub(crate) fn ramachandran(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyRamachandranOptions,
    outliers_only: bool,
) -> PyResult<Vec<PyRamachandranRecord>> {
    project_ramachandran(py, structure, options, outliers_only)
}

#[pyfunction]
pub(crate) fn ramachandran_outliers(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyRamachandranOptions,
) -> PyResult<Vec<PyRamachandranRecord>> {
    project_ramachandran(py, structure, options, true)
}

fn project_ramachandran(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyRamachandranOptions,
    outliers_only: bool,
) -> PyResult<Vec<PyRamachandranRecord>> {
    let structure = structure.structure().clone();
    let references = options.references.clone();
    let basins = options.basins.clone();
    let minimum_probability = options.minimum_probability;
    py.detach(move || {
        let options =
            pdbiox::validate::RamachandranOptions::new(&references, basins, minimum_probability)
                .map_err(|error| error.to_string())?;
        if outliers_only {
            pdbiox::validate::ramachandran_outliers(&structure, &options)
        } else {
            pdbiox::validate::ramachandran(&structure, &options)
        }
        .map_err(|error| error.to_string())
    })
    .map(|values| values.into_iter().map(PyRamachandranRecord::from).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn altloc_occupancy_sums(
    py: Python<'_>,
    structure: &PyStructure,
    namespace: PyNamespace,
    options: PyAltlocOccupancyOptions,
) -> PyResult<super::PyAltlocOccupancyReport> {
    let structure = structure.structure().clone();
    py.detach(move || {
        pdbiox::validate::altloc_occupancy_sums(&structure, namespace.into(), options.0)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn ccd_missing_atoms(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &RootComponentDictionary,
    policy: &PyAnalysisPolicy,
) -> PyResult<super::PyCcdCompletenessReport> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let policy = policy.inner.clone();
    py.detach(move || pdbiox::validate::ccd_missing_atoms(&structure, dictionary.as_ref(), &policy))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn plane_restraint_outliers(
    py: Python<'_>,
    structure: &PyStructure,
    restraints: Vec<super::PyPlaneRestraint>,
    options: PyPlanarityOptions,
) -> PyResult<super::PyPlaneRestraintReport> {
    let structure = structure.structure().clone();
    let restraints = restraints
        .into_iter()
        .map(|value| pdbiox::validate::PlaneRestraint {
            id: value.id,
            atoms: value.atoms.inner,
        })
        .collect::<Vec<_>>();
    py.detach(move || {
        pdbiox::validate::plane_restraint_outliers(&structure, &restraints, options.0)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, dictionary, options, *, policy=None))]
pub(crate) fn chirality_outliers(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &RootComponentDictionary,
    options: &PyChiralityOptions,
    policy: Option<&PyAnalysisPolicy>,
) -> PyResult<PyChiralityReport> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let options = options.0;
    let policy = policy.map_or_else(pdbiox::AnalysisPolicy::default, |value| value.inner.clone());
    py.detach(move || {
        pdbiox::validate::chirality_outliers(&structure, dictionary.as_ref(), &policy, options)
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, dictionary, references, profile, options, *, policy=None))]
pub(crate) fn rotamer_outliers(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &RootComponentDictionary,
    references: &PyReferenceLibrary,
    profile: &PyRotamerProfile,
    options: &PyRotamerOptions,
    policy: Option<&PyAnalysisPolicy>,
) -> PyResult<PyRotamerReport> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let references = references.0.clone();
    let profile = profile.0.clone();
    let options = options.0;
    let policy = policy.map_or_else(pdbiox::AnalysisPolicy::default, |value| value.inner.clone());
    py.detach(move || {
        pdbiox::validate::rotamer_outliers(
            &structure,
            dictionary.as_ref(),
            &policy,
            &references,
            &profile,
            options,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn classify(
    py: Python<'_>,
    phi: f64,
    psi: f64,
    options: &PyRamachandranOptions,
) -> PyResult<(super::PyRamachandranRegion, super::PyReferenceAssessment)> {
    let basins = options.basins.clone();
    let references = options.references.clone();
    let minimum_probability = options.minimum_probability;
    py.detach(move || {
        let options =
            pdbiox::validate::RamachandranOptions::new(&references, basins, minimum_probability)?;
        pdbiox::validate::classify(phi, psi, &options)
    })
    .map(|(region, assessment)| (region.into(), assessment.into()))
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn reference_geometry(
    py: Python<'_>,
    structure: &PyStructure,
    provider: &RootComponentDictionary,
    namespace: PyNamespace,
    options: PyReferenceGeometryOptions,
) -> PyResult<PyReferenceGeometryReport> {
    let structure = structure.structure().clone();
    let provider = provider.0.clone();
    py.detach(move || {
        pdbiox::validate::reference_geometry(
            &structure,
            provider.as_ref(),
            namespace.into(),
            options.0,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn nucleic_acid_geometry(
    py: Python<'_>,
    structure: &PyStructure,
    provider: &RootComponentDictionary,
    roles: &PyPolymerRoleProfile,
    policy: PyNucleicGeometryPolicy,
) -> PyResult<Vec<PyNucleicGeometryRecord>> {
    let structure = structure.structure().clone();
    let provider = provider.0.clone();
    let roles = roles.0.clone();
    py.detach(move || {
        pdbiox::validate::nucleic_acid_geometry(&structure, provider.as_ref(), &roles, policy.0)
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(value_error)
}
