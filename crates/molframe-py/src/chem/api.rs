//! Free chemistry operations that accept the same native provider primitives.

use super::components::{PyComponent, PyComponentCoverage};
use super::providers::extract_provider;
use super::roles::{
    PyChemistryReport, PyPolymerLinkPolicy, PyPolymerRoleProfile, PyPolymerRoleReport,
};
use super::{PeoeError, PyPeoeOptions};
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use ::molframe;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn component_coverage(
    py: Python<'_>,
    structure: &PyStructure,
    provider: &Bound<'_, PyAny>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyComponentCoverage> {
    let structure = structure.structure().clone();
    let provider = extract_provider(provider)?;
    let policy = policy.inner.clone();
    py.detach(move || molframe::component_coverage(&structure, &provider, &policy))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn component_peoe_charges(
    py: Python<'_>,
    component: &PyComponent,
    options: PyPeoeOptions,
) -> PyResult<Vec<f64>> {
    let component = component.0.clone();
    py.detach(move || molframe::component_peoe_charges(&component, options.0))
        .map_err(|error| PeoeError::new_err(error.to_string()))
}

#[pyfunction]
#[pyo3(signature = (structure, provider, policy=None))]
pub(crate) fn apply_component_chemistry(
    py: Python<'_>,
    structure: &PyStructure,
    provider: &Bound<'_, PyAny>,
    policy: Option<PyPolymerLinkPolicy>,
) -> PyResult<PyChemistryReport> {
    let structure = structure.structure().clone();
    let provider = extract_provider(provider)?;
    let policy = policy.map_or(molframe::PolymerLinkPolicy::Disabled, |value| value.0);
    py.detach(move || molframe::apply_component_chemistry(&structure, &provider, policy))
        .map(Into::into)
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn apply_polymer_role_profile(
    py: Python<'_>,
    structure: &PyStructure,
    provider: &Bound<'_, PyAny>,
    profile: &PyPolymerRoleProfile,
) -> PyResult<PyPolymerRoleReport> {
    let structure = structure.structure().clone();
    let provider = extract_provider(provider)?;
    let profile = profile.0.clone();
    py.detach(move || molframe::apply_polymer_role_profile(&structure, &provider, &profile))
        .map(Into::into)
        .map_err(value_error)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
