//! Python projections for strict application and analysis-policy configuration.

use crate::query::PyAnalysisPolicy;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyclass(name = "PolicyOverrides", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolicyOverrides {
    #[pyo3(get)]
    assembly: Option<String>,
    #[pyo3(get)]
    model: Option<String>,
    #[pyo3(get)]
    altloc: Option<String>,
    #[pyo3(get)]
    identifiers: Option<String>,
    #[pyo3(get)]
    missing_atoms: Option<String>,
    #[pyo3(get)]
    hydrogens: Option<String>,
    #[pyo3(get)]
    atom_equivalence: Option<String>,
    #[pyo3(get)]
    symmetry: Option<String>,
    #[pyo3(get)]
    alignment: Option<String>,
    #[pyo3(get)]
    precision: Option<String>,
    #[pyo3(get)]
    periodic: Option<String>,
    #[pyo3(get)]
    vdw_radii: Option<String>,
    #[pyo3(get)]
    contact_def: Option<String>,
    #[pyo3(get)]
    float_tolerance_relative: Option<f64>,
    #[pyo3(get)]
    float_tolerance_absolute: Option<f64>,
}

#[pyclass(name = "OutputConfiguration", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyOutputConfiguration {
    #[pyo3(get)]
    format: Option<String>,
}

#[pyclass(name = "ChemistryConfiguration", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChemistryConfiguration {
    #[pyo3(get)]
    ccd_cache: Option<PathBuf>,
    #[pyo3(get)]
    ccd_version: Option<String>,
}

#[pyclass(name = "ApplicationConfiguration", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyApplicationConfiguration {
    #[pyo3(get)]
    policy: PyPolicyOverrides,
    #[pyo3(get)]
    output: PyOutputConfiguration,
    #[pyo3(get)]
    chem: PyChemistryConfiguration,
}

#[pymethods]
impl PyPolicyOverrides {
    #[new]
    #[expect(
        clippy::too_many_arguments,
        reason = "The Python constructor intentionally mirrors the public policy field set."
    )]
    #[pyo3(signature = (*, assembly=None, model=None, altloc=None, identifiers=None,
        missing_atoms=None, hydrogens=None, atom_equivalence=None, symmetry=None,
        alignment=None, precision=None, periodic=None, vdw_radii=None,
        contact_def=None, float_tolerance_relative=None, float_tolerance_absolute=None))]
    fn new(
        assembly: Option<String>,
        model: Option<String>,
        altloc: Option<String>,
        identifiers: Option<String>,
        missing_atoms: Option<String>,
        hydrogens: Option<String>,
        atom_equivalence: Option<String>,
        symmetry: Option<String>,
        alignment: Option<String>,
        precision: Option<String>,
        periodic: Option<String>,
        vdw_radii: Option<String>,
        contact_def: Option<String>,
        float_tolerance_relative: Option<f64>,
        float_tolerance_absolute: Option<f64>,
    ) -> Self {
        Self {
            assembly,
            model,
            altloc,
            identifiers,
            missing_atoms,
            hydrogens,
            atom_equivalence,
            symmetry,
            alignment,
            precision,
            periodic,
            vdw_radii,
            contact_def,
            float_tolerance_relative,
            float_tolerance_absolute,
        }
    }

    #[pyo3(signature = (policy=None))]
    fn apply_to(&self, policy: Option<&PyAnalysisPolicy>) -> PyResult<PyAnalysisPolicy> {
        let overrides = self.clone().into_inner();
        let base = policy.map_or_else(molframe::AnalysisPolicy::default, |value| {
            value.inner.clone()
        });
        overrides
            .apply_to(base)
            .map(|inner| PyAnalysisPolicy { inner })
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}

#[pymethods]
impl PyOutputConfiguration {
    #[new]
    #[pyo3(signature = (format=None))]
    fn new(format: Option<String>) -> Self {
        Self { format }
    }
}

#[pymethods]
impl PyChemistryConfiguration {
    #[new]
    #[pyo3(signature = (ccd_cache=None, ccd_version=None))]
    fn new(ccd_cache: Option<PathBuf>, ccd_version: Option<String>) -> Self {
        Self {
            ccd_cache,
            ccd_version,
        }
    }
}

#[pyfunction]
pub(crate) fn read_configuration(
    py: Python<'_>,
    path: PathBuf,
) -> PyResult<PyApplicationConfiguration> {
    py.detach(move || molframe::read_configuration(path))
        .map(Into::into)
        .map_err(|error| crate::errors::PolicyConfigError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn read_policy_overrides(py: Python<'_>, path: PathBuf) -> PyResult<PyPolicyOverrides> {
    py.detach(move || molframe::read_policy_overrides(path))
        .map(Into::into)
        .map_err(|error| crate::errors::PolicyConfigError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn read_policy(py: Python<'_>, path: PathBuf) -> PyResult<PyAnalysisPolicy> {
    py.detach(move || molframe::read_policy(path))
        .map(|inner| PyAnalysisPolicy { inner })
        .map_err(|error| crate::errors::PolicyConfigError::new_err(error.to_string()))
}

impl PyPolicyOverrides {
    fn into_inner(self) -> molframe::PolicyOverrides {
        molframe::PolicyOverrides {
            assembly: self.assembly,
            model: self.model,
            altloc: self.altloc,
            identifiers: self.identifiers,
            missing_atoms: self.missing_atoms,
            hydrogens: self.hydrogens,
            atom_equivalence: self.atom_equivalence,
            symmetry: self.symmetry,
            alignment: self.alignment,
            precision: self.precision,
            periodic: self.periodic,
            vdw_radii: self.vdw_radii,
            contact_def: self.contact_def,
            float_tolerance_relative: self.float_tolerance_relative,
            float_tolerance_absolute: self.float_tolerance_absolute,
        }
    }
}

impl From<molframe::PolicyOverrides> for PyPolicyOverrides {
    fn from(value: molframe::PolicyOverrides) -> Self {
        Self {
            assembly: value.assembly,
            model: value.model,
            altloc: value.altloc,
            identifiers: value.identifiers,
            missing_atoms: value.missing_atoms,
            hydrogens: value.hydrogens,
            atom_equivalence: value.atom_equivalence,
            symmetry: value.symmetry,
            alignment: value.alignment,
            precision: value.precision,
            periodic: value.periodic,
            vdw_radii: value.vdw_radii,
            contact_def: value.contact_def,
            float_tolerance_relative: value.float_tolerance_relative,
            float_tolerance_absolute: value.float_tolerance_absolute,
        }
    }
}

impl From<molframe::ApplicationConfiguration> for PyApplicationConfiguration {
    fn from(value: molframe::ApplicationConfiguration) -> Self {
        Self {
            policy: value.policy.into(),
            output: value.output.into(),
            chem: value.chem.into(),
        }
    }
}

impl From<molframe::OutputConfiguration> for PyOutputConfiguration {
    fn from(value: molframe::OutputConfiguration) -> Self {
        Self {
            format: value.format,
        }
    }
}

impl From<molframe::ChemistryConfiguration> for PyChemistryConfiguration {
    fn from(value: molframe::ChemistryConfiguration) -> Self {
        Self {
            ccd_cache: value.ccd_cache,
            ccd_version: value.ccd_version,
        }
    }
}

#[cfg(test)]
#[path = "config_tests.rs"]
mod tests;
