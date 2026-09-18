//! Python-owned policies for applying CCD chemistry and polymer roles.

use super::PyComponentDictionary;
use super::components::{PyComponentKind, PyPolymerAtomRole};
use crate::contract::PyDiagnostic;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "PolymerLinkRule", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolymerLinkRule {
    #[pyo3(get)]
    left_kind: PyComponentKind,
    #[pyo3(get)]
    left_atom: String,
    #[pyo3(get)]
    right_kind: PyComponentKind,
    #[pyo3(get)]
    right_atom: String,
}

#[pymethods]
impl PyPolymerLinkRule {
    #[new]
    fn new(
        left_kind: PyComponentKind,
        left_atom: String,
        right_kind: PyComponentKind,
        right_atom: String,
    ) -> Self {
        Self {
            left_kind,
            left_atom,
            right_kind,
            right_atom,
        }
    }
}

#[pyclass(name = "PolymerLinkPolicy", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolymerLinkPolicy(pub(crate) molframe::PolymerLinkPolicy);

#[pymethods]
impl PyPolymerLinkPolicy {
    #[staticmethod]
    fn disabled() -> Self {
        Self(molframe::PolymerLinkPolicy::Disabled)
    }

    #[staticmethod]
    fn explicit(angstrom: f32, rules: Vec<PyPolymerLinkRule>) -> Self {
        let rules = rules.into_iter().map(Into::into).collect::<Vec<_>>();
        Self(molframe::PolymerLinkPolicy::explicit(angstrom, rules))
    }

    #[getter]
    fn is_disabled(&self) -> bool {
        matches!(self.0, molframe::PolymerLinkPolicy::Disabled)
    }
}

#[pyclass(name = "ChemistryReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChemistryReport {
    #[pyo3(get)]
    structure: PyStructure,
    #[pyo3(get)]
    findings: Vec<PyDiagnostic>,
    #[pyo3(get)]
    dictionary_version: String,
    #[pyo3(get)]
    polymer_link_policy: PyPolymerLinkPolicy,
}

#[pyclass(name = "PolymerRoleRule", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolymerRoleRule {
    #[pyo3(get)]
    component_id: Option<String>,
    #[pyo3(get)]
    component_kind: Option<PyComponentKind>,
    #[pyo3(get)]
    atom_name: String,
    #[pyo3(get)]
    role: PyPolymerAtomRole,
}

#[pymethods]
impl PyPolymerRoleRule {
    #[new]
    fn new(
        component_id: Option<String>,
        component_kind: Option<PyComponentKind>,
        atom_name: String,
        role: PyPolymerAtomRole,
    ) -> Self {
        Self {
            component_id,
            component_kind,
            atom_name,
            role,
        }
    }
}

#[pyclass(name = "PolymerRoleProfile", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolymerRoleProfile(pub(crate) molframe::PolymerRoleProfile);

#[pymethods]
impl PyPolymerRoleProfile {
    #[new]
    fn new(id: String, rules: Vec<PyPolymerRoleRule>) -> Self {
        Self(molframe::PolymerRoleProfile {
            id: id.into(),
            rules: rules.into_iter().map(Into::into).collect::<Vec<_>>().into(),
        })
    }

    #[getter]
    fn id(&self) -> String {
        self.0.id.to_string()
    }

    #[getter]
    fn rules(&self) -> Vec<PyPolymerRoleRule> {
        self.0.rules.iter().cloned().map(Into::into).collect()
    }
}

#[pyclass(name = "PolymerRoleReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolymerRoleReport {
    #[pyo3(get)]
    structure: PyStructure,
    #[pyo3(get)]
    unresolved_components: Vec<String>,
    #[pyo3(get)]
    dictionary_version: String,
    #[pyo3(get)]
    profile_id: String,
}

#[pymethods]
impl PyComponentDictionary {
    #[pyo3(signature = (structure, policy=None))]
    fn apply_chemistry(
        &self,
        structure: &PyStructure,
        policy: Option<PyPolymerLinkPolicy>,
    ) -> PyResult<PyChemistryReport> {
        molframe::apply_component_chemistry(
            structure.structure(),
            self.0.as_ref(),
            policy.map_or(molframe::PolymerLinkPolicy::Disabled, |value| value.0),
        )
        .map(Into::into)
        .map_err(value_error)
    }

    fn apply_polymer_role_profile(
        &self,
        structure: &PyStructure,
        profile: &PyPolymerRoleProfile,
    ) -> PyResult<PyPolymerRoleReport> {
        molframe::apply_polymer_role_profile(structure.structure(), self.0.as_ref(), &profile.0)
            .map(Into::into)
            .map_err(value_error)
    }

    fn coverage_with_policy(
        &self,
        structure: &PyStructure,
        policy: &PyAnalysisPolicy,
    ) -> PyResult<super::components::PyComponentCoverage> {
        self.coverage(structure, policy)
    }
}

impl From<PyPolymerLinkRule> for molframe::PolymerLinkRule {
    fn from(value: PyPolymerLinkRule) -> Self {
        Self::new(
            value.left_kind.into(),
            value.left_atom,
            value.right_kind.into(),
            value.right_atom,
        )
    }
}

impl From<PyPolymerRoleRule> for molframe::PolymerRoleRule {
    fn from(value: PyPolymerRoleRule) -> Self {
        Self {
            component_id: value.component_id.map(Into::into),
            component_kind: value.component_kind.map(Into::into),
            atom_name: value.atom_name.into(),
            role: value.role.0,
        }
    }
}

impl From<molframe::PolymerLinkPolicy> for PyPolymerLinkPolicy {
    fn from(value: molframe::PolymerLinkPolicy) -> Self {
        Self(value)
    }
}

impl From<molframe::ChemistryReport> for PyChemistryReport {
    fn from(value: molframe::ChemistryReport) -> Self {
        Self {
            structure: PyStructure::new(value.structure),
            findings: value.findings.into_iter().map(Into::into).collect(),
            dictionary_version: value.dictionary_version.as_str().to_owned(),
            polymer_link_policy: value.polymer_link_policy.into(),
        }
    }
}

impl From<molframe::PolymerRoleRule> for PyPolymerRoleRule {
    fn from(value: molframe::PolymerRoleRule) -> Self {
        Self {
            component_id: value.component_id.map(|value| value.to_string()),
            component_kind: value.component_kind.map(Into::into),
            atom_name: value.atom_name.to_string(),
            role: PyPolymerAtomRole(value.role),
        }
    }
}

impl From<molframe::PolymerRoleReport> for PyPolymerRoleReport {
    fn from(value: molframe::PolymerRoleReport) -> Self {
        Self {
            structure: PyStructure::new(value.structure),
            unresolved_components: value
                .unresolved_components
                .into_iter()
                .map(|value| value.to_string())
                .collect(),
            dictionary_version: value.dictionary_version.as_str().to_owned(),
            profile_id: value.profile_id.to_string(),
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPolymerLinkRule>()?;
    module.add_class::<PyPolymerLinkPolicy>()?;
    module.add_class::<PyChemistryReport>()?;
    module.add_class::<PyPolymerRoleRule>()?;
    module.add_class::<PyPolymerRoleProfile>()?;
    module.add_class::<PyPolymerRoleReport>()?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
