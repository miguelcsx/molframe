//! Conversion boundary for caller-authored CCD polymer role profiles.

use crate::bindings::{PyStructure, findings_error};
use molframe::chemistry as chem;
use pyo3::prelude::*;
use std::path::PathBuf;

#[derive(Clone)]
#[pyclass(
    name = "PolymerRoleRule",
    frozen,
    skip_from_py_object,
    module = "molframe.chemistry"
)]
struct PyPolymerRoleRule {
    inner: chem::PolymerRoleRule,
}

#[pymethods]
impl PyPolymerRoleRule {
    #[new]
    #[pyo3(signature = (atom_name, role, *, component_id=None, component_kind=None))]
    fn new(
        atom_name: String,
        role: i64,
        component_id: Option<String>,
        component_kind: Option<i64>,
    ) -> PyResult<Self> {
        let role = chem::PolymerAtomRole::from_code(role)
            .ok_or_else(|| crate::error::value("invalid polymer role bits"))?;
        let component_kind = component_kind
            .map(|code| {
                chem::ComponentKind::from_code(code)
                    .ok_or_else(|| crate::error::value("invalid component kind code"))
            })
            .transpose()?;
        Ok(Self {
            inner: chem::PolymerRoleRule {
                atom_name: atom_name.into(),
                role,
                component_id: component_id.map(String::into_boxed_str),
                component_kind,
            },
        })
    }

    #[getter]
    fn atom_name(&self) -> &str {
        &self.inner.atom_name
    }
    #[getter]
    fn role(&self) -> i64 {
        self.inner.role.code()
    }
    #[getter]
    fn component_id(&self) -> Option<&str> {
        self.inner.component_id.as_deref()
    }
    #[getter]
    fn component_kind(&self) -> Option<i64> {
        self.inner.component_kind.map(chem::ComponentKind::code)
    }
}

#[pyclass(
    name = "PolymerRoleReport",
    frozen,
    skip_from_py_object,
    module = "molframe.chemistry"
)]
struct PyPolymerRoleReport {
    #[pyo3(get)]
    structure: PyStructure,
    #[pyo3(get)]
    unresolved_components: Vec<String>,
    #[pyo3(get)]
    dictionary_version: String,
    #[pyo3(get)]
    profile_id: String,
}

#[pyfunction]
#[pyo3(signature = (structure, components, rules, *, profile_id, version="unversioned"))]
fn apply_polymer_role_profile(
    py: Python<'_>,
    structure: &PyStructure,
    components: PathBuf,
    rules: Vec<PyRef<'_, PyPolymerRoleRule>>,
    profile_id: String,
    version: &str,
) -> PyResult<PyPolymerRoleReport> {
    let source = structure.inner.clone();
    let profile = chem::PolymerRoleProfile {
        id: profile_id.into(),
        rules: rules.into_iter().map(|rule| rule.inner.clone()).collect(),
    };
    let version = molframe::DictionaryVersion::new(version);
    let report = py.detach(move || {
        let (provider, _) = molframe::read_component_dictionary(&components, version)
            .map_err(|findings| findings_error(&findings))?;
        chem::apply_polymer_role_profile(source.engine(), &provider, &profile)
            .map_err(crate::error::kernel)
    })?;
    Ok(PyPolymerRoleReport {
        structure: PyStructure::new(molframe::Structure::from(report.structure)),
        unresolved_components: report
            .unresolved_components
            .into_iter()
            .map(String::from)
            .collect(),
        dictionary_version: report.dictionary_version.as_str().to_owned(),
        profile_id: report.profile_id.into_string(),
    })
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPolymerRoleRule>()?;
    module.add_class::<PyPolymerRoleReport>()?;
    module.add_function(wrap_pyfunction!(apply_polymer_role_profile, module)?)
}
