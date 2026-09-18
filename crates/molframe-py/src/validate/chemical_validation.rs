//! CCD- and reference-governed stereochemistry validation bindings.

use crate::analysis::PyReferenceLibrary;
use crate::chemistry::PyComponentDictionary;
use crate::chemistry::components::PyStereoConfiguration;
use crate::contract::PyDiagnostic;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "ChiralityIssue", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyChiralityIssue {
    Inverted,
    Degenerate,
}

#[pyclass(name = "ChiralityOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChiralityOptions(pub(crate) molframe::validate::ChiralityOptions);

#[pymethods]
impl PyChiralityOptions {
    #[new]
    fn new(minimum_abs_volume: f64) -> Self {
        Self(molframe::validate::ChiralityOptions { minimum_abs_volume })
    }
}

#[pyclass(name = "ChiralityFlag", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChiralityFlag {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    centre: u32,
    #[pyo3(get)]
    expected: PyStereoConfiguration,
    #[pyo3(get)]
    issue: PyChiralityIssue,
    #[pyo3(get)]
    observed_volume: f64,
    #[pyo3(get)]
    reference_volume: f64,
}

#[pyclass(name = "ChiralityReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChiralityReport {
    #[pyo3(get)]
    flags: Vec<PyChiralityFlag>,
    #[pyo3(get)]
    findings: Vec<PyDiagnostic>,
    #[pyo3(get)]
    dictionary_version: String,
}

#[pyclass(name = "RotamerDefinition", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRotamerDefinition(molframe::validate::RotamerDefinition);

#[pymethods]
impl PyRotamerDefinition {
    #[new]
    fn new(component: &str, chi_index: u8, atoms: [String; 4], distribution: &str) -> Self {
        Self(molframe::validate::RotamerDefinition {
            component_id: component.into(),
            chi_index,
            atoms: atoms.map(Into::into),
            distribution: distribution.into(),
        })
    }
}

#[pyclass(name = "RotamerProfile", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRotamerProfile(pub(super) molframe::validate::RotamerProfile);

#[pymethods]
impl PyRotamerProfile {
    #[new]
    fn new(id: &str, version: &str, definitions: Vec<PyRotamerDefinition>) -> PyResult<Self> {
        molframe::validate::RotamerProfile::new(
            id,
            version,
            definitions.into_iter().map(|value| value.0),
        )
        .map(Self)
        .map_err(value_error)
    }

    #[getter]
    fn id(&self) -> &str {
        self.0.id()
    }
    #[getter]
    fn version(&self) -> &str {
        self.0.version()
    }
}

#[pyclass(name = "RotamerOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRotamerOptions(pub(crate) molframe::validate::RotamerOptions);

#[pymethods]
impl PyRotamerOptions {
    #[new]
    fn new(minimum_probability: f64) -> Self {
        Self(molframe::validate::RotamerOptions {
            minimum_probability,
        })
    }
}

#[pyclass(name = "RotamerFlag", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRotamerFlag {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    chi_index: u8,
    #[pyo3(get)]
    chi_degrees: f64,
    #[pyo3(get)]
    probability: f64,
    #[pyo3(get)]
    distribution: String,
}

#[pyclass(name = "RotamerReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRotamerReport {
    #[pyo3(get)]
    flags: Vec<PyRotamerFlag>,
    #[pyo3(get)]
    findings: Vec<PyDiagnostic>,
    #[pyo3(get)]
    dictionary_version: String,
    #[pyo3(get)]
    profile_id: String,
    #[pyo3(get)]
    profile_version: String,
    #[pyo3(get)]
    reference_id: String,
    #[pyo3(get)]
    reference_version: String,
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (dictionary, options, *, policy=None))]
    fn chirality_outliers(
        &self,
        py: Python<'_>,
        dictionary: &PyComponentDictionary,
        options: &PyChiralityOptions,
        policy: Option<&PyAnalysisPolicy>,
    ) -> PyResult<PyChiralityReport> {
        let structure = self.structure().clone();
        let dictionary = dictionary.0.clone();
        let options = options.0;
        let policy = policy.map_or_else(molframe::AnalysisPolicy::default, |value| {
            value.inner.clone()
        });
        py.detach(move || {
            molframe::validate::chirality_outliers(
                &structure,
                dictionary.as_ref(),
                &policy,
                options,
            )
        })
        .map(PyChiralityReport::from)
        .map_err(value_error)
    }

    #[pyo3(signature = (dictionary, references, profile, options, *, policy=None))]
    fn rotamer_outliers(
        &self,
        py: Python<'_>,
        dictionary: &PyComponentDictionary,
        references: &PyReferenceLibrary,
        profile: &PyRotamerProfile,
        options: &PyRotamerOptions,
        policy: Option<&PyAnalysisPolicy>,
    ) -> PyResult<PyRotamerReport> {
        let structure = self.structure().clone();
        let dictionary = dictionary.0.clone();
        let references = references.0.clone();
        let profile = profile.0.clone();
        let options = options.0;
        let policy = policy.map_or_else(molframe::AnalysisPolicy::default, |value| {
            value.inner.clone()
        });
        py.detach(move || {
            molframe::validate::rotamer_outliers(
                &structure,
                dictionary.as_ref(),
                &policy,
                &references,
                &profile,
                options,
            )
        })
        .map(PyRotamerReport::from)
        .map_err(value_error)
    }
}

impl From<molframe::validate::ChiralityReport> for PyChiralityReport {
    fn from(value: molframe::validate::ChiralityReport) -> Self {
        Self {
            flags: value.flags.into_iter().map(PyChiralityFlag::from).collect(),
            findings: value.findings.into_iter().map(Into::into).collect(),
            dictionary_version: value.dictionary_version.as_str().to_owned(),
        }
    }
}

impl From<molframe::validate::ChiralityFlag> for PyChiralityFlag {
    fn from(value: molframe::validate::ChiralityFlag) -> Self {
        Self {
            residue: value.residue.get(),
            centre: value.centre.get(),
            expected: value.expected.into(),
            issue: value.issue.into(),
            observed_volume: value.observed_volume,
            reference_volume: value.reference_volume,
        }
    }
}

impl From<molframe::validate::ChiralityIssue> for PyChiralityIssue {
    fn from(value: molframe::validate::ChiralityIssue) -> Self {
        match value {
            molframe::validate::ChiralityIssue::Inverted => Self::Inverted,
            molframe::validate::ChiralityIssue::Degenerate => Self::Degenerate,
        }
    }
}

impl From<molframe::validate::RotamerReport> for PyRotamerReport {
    fn from(value: molframe::validate::RotamerReport) -> Self {
        Self {
            flags: value.flags.into_iter().map(PyRotamerFlag::from).collect(),
            findings: value.findings.into_iter().map(Into::into).collect(),
            dictionary_version: value.dictionary_version.as_str().to_owned(),
            profile_id: value.profile_id.into(),
            profile_version: value.profile_version.into(),
            reference_id: value.reference_id.into(),
            reference_version: value.reference_version.into(),
        }
    }
}

impl From<molframe::validate::RotamerFlag> for PyRotamerFlag {
    fn from(value: molframe::validate::RotamerFlag) -> Self {
        Self {
            residue: value.residue.get(),
            chi_index: value.chi_index,
            chi_degrees: value.chi_degrees,
            probability: value.probability,
            distribution: value.distribution.into(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
