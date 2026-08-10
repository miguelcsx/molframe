//! Typed controls for the policy dimensions used by selection.

use pdbiox::{AltlocPolicy, AnalysisPolicy, AssemblyChoice, MissingPolicy, ModelChoice, Namespace};
use pyo3::prelude::*;

macro_rules! policy_enum {
    ($python:literal, $name:ident { $($variant:ident),+ $(,)? }) => {
        #[pyclass(name = $python, frozen, eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) enum $name { $($variant),+ }
    };
}

policy_enum!(
    "Namespace",
    PyNamespace {
        Label,
        Auth,
        Explicit
    }
);
policy_enum!(
    "MissingPolicy",
    PyMissingPolicy {
        Ignore,
        Report,
        Indeterminate,
        Fail
    }
);

#[pyclass(name = "ModelChoice", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyModelChoice(ModelChoice);

#[pymethods]
impl PyModelChoice {
    #[staticmethod]
    fn first() -> Self {
        Self(ModelChoice::First)
    }
    #[staticmethod]
    fn index(index: u32) -> Self {
        Self(ModelChoice::Index(index))
    }
    #[staticmethod]
    fn all() -> Self {
        Self(ModelChoice::All)
    }
    #[staticmethod]
    fn ensemble() -> Self {
        Self(ModelChoice::Ensemble)
    }
}

#[pyclass(name = "AltlocPolicy", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAltlocPolicy(AltlocPolicy);

#[pymethods]
impl PyAltlocPolicy {
    #[staticmethod]
    fn keep_all() -> Self {
        Self(AltlocPolicy::KeepAll)
    }
    #[staticmethod]
    fn conformer_consistent() -> Self {
        Self(AltlocPolicy::ConformerConsistent)
    }
    #[staticmethod]
    fn label(label: String) -> Self {
        Self(AltlocPolicy::Label(label.into()))
    }
    #[staticmethod]
    fn first() -> Self {
        Self(AltlocPolicy::First)
    }
    #[staticmethod]
    fn highest_occupancy_per_residue() -> Self {
        Self(AltlocPolicy::HighestOccupancyPerResidue)
    }
    #[staticmethod]
    fn highest_occupancy_per_atom() -> Self {
        Self(AltlocPolicy::HighestOccupancyPerAtom)
    }
}

#[pyclass(name = "AssemblyChoice", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAssemblyChoice(AssemblyChoice);

#[pymethods]
impl PyAssemblyChoice {
    #[staticmethod]
    fn asymmetric_unit() -> Self {
        Self(AssemblyChoice::AsymmetricUnit)
    }
    #[staticmethod]
    fn biological(identifier: String) -> Self {
        Self(AssemblyChoice::Biological(identifier.into()))
    }
    #[staticmethod]
    fn crystal(radius: f32) -> Self {
        Self(AssemblyChoice::Crystal { radius })
    }
}

#[pyclass(name = "AnalysisPolicy", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAnalysisPolicy {
    pub(crate) inner: AnalysisPolicy,
}

#[pymethods]
impl PyAnalysisPolicy {
    #[new]
    #[pyo3(signature = (*, assembly=None, model=None, altloc=None, identifiers=PyNamespace::Auth, missing_atoms=PyMissingPolicy::Report))]
    fn new(
        assembly: Option<&PyAssemblyChoice>,
        model: Option<&PyModelChoice>,
        altloc: Option<&PyAltlocPolicy>,
        identifiers: PyNamespace,
        missing_atoms: PyMissingPolicy,
    ) -> Self {
        let mut inner = AnalysisPolicy::default()
            .with_identifiers(identifiers.into())
            .with_missing_atoms(missing_atoms.into());
        if let Some(value) = assembly {
            inner = inner.with_assembly(value.0.clone());
        }
        if let Some(value) = model {
            inner = inner.with_model(value.0);
        }
        if let Some(value) = altloc {
            inner = inner.with_altloc(value.0.clone());
        }
        Self { inner }
    }

    #[getter]
    fn fingerprint(&self) -> String {
        self.inner.fingerprint().to_string()
    }
}

impl From<PyNamespace> for Namespace {
    fn from(value: PyNamespace) -> Self {
        match value {
            PyNamespace::Label => Self::Label,
            PyNamespace::Auth => Self::Auth,
            PyNamespace::Explicit => Self::Explicit,
        }
    }
}

impl From<PyMissingPolicy> for MissingPolicy {
    fn from(value: PyMissingPolicy) -> Self {
        match value {
            PyMissingPolicy::Ignore => Self::Ignore,
            PyMissingPolicy::Report => Self::Report,
            PyMissingPolicy::Indeterminate => Self::Indeterminate,
            PyMissingPolicy::Fail => Self::Fail,
        }
    }
}
