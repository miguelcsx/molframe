//! Typed controls for the policy dimensions used by selection.

use crate::audit::{
    PyAlignmentPolicy, PyContactDefinition, PyEquivalencePolicy, PyHydrogenPolicy,
    PyPeriodicPolicy, PyPrecision, PySymmetryPolicy, PyTolerance,
};
use crate::chemistry::PyRadiusSet;
use molframe::{
    AltlocPolicy, AnalysisPolicy, AssemblyChoice, MissingPolicy, ModelChoice, Namespace,
};
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
pub(crate) struct PyModelChoice(pub(crate) ModelChoice);

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
pub(crate) struct PyAltlocPolicy(pub(crate) AltlocPolicy);

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
pub(crate) struct PyAssemblyChoice(pub(crate) AssemblyChoice);

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
    #[pyo3(signature = (*, assembly=None, model=None, altloc=None,
        identifiers=PyNamespace::Auth, missing_atoms=PyMissingPolicy::Report,
        hydrogens=None, atom_equivalence=None, symmetry=None, alignment=None,
        precision=None, periodic=None, vdw_radii=None, contact_def=None,
        float_tolerance=None))]
    fn new(
        assembly: Option<&PyAssemblyChoice>,
        model: Option<&PyModelChoice>,
        altloc: Option<&PyAltlocPolicy>,
        identifiers: PyNamespace,
        missing_atoms: PyMissingPolicy,
        hydrogens: Option<PyHydrogenPolicy>,
        atom_equivalence: Option<PyEquivalencePolicy>,
        symmetry: Option<PySymmetryPolicy>,
        alignment: Option<PyAlignmentPolicy>,
        precision: Option<PyPrecision>,
        periodic: Option<PyPeriodicPolicy>,
        vdw_radii: Option<PyRadiusSet>,
        contact_def: Option<PyContactDefinition>,
        float_tolerance: Option<PyTolerance>,
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
        if let Some(value) = hydrogens {
            inner.hydrogens = value.into();
        }
        if let Some(value) = atom_equivalence {
            inner.atom_equivalence = value.into();
        }
        if let Some(value) = symmetry {
            inner.symmetry = value.into();
        }
        if let Some(value) = alignment {
            inner.alignment = value.into();
        }
        if let Some(value) = precision {
            inner.precision = value.into();
        }
        if let Some(value) = periodic {
            inner.periodic = value.into();
        }
        if let Some(value) = vdw_radii {
            inner.vdw_radii = match value {
                PyRadiusSet::Bondi => molframe::core::contract::RadiiSet::Bondi,
                PyRadiusSet::AmberUnited => molframe::core::contract::RadiiSet::AmberUnited,
                PyRadiusSet::Charmm => molframe::core::contract::RadiiSet::Charmm,
                PyRadiusSet::Alvarez => molframe::core::contract::RadiiSet::Alvarez,
            };
        }
        if let Some(value) = contact_def {
            inner.contact_def = value.into();
        }
        if let Some(value) = float_tolerance {
            inner.float_tolerance = value.into();
        }
        Self { inner }
    }

    #[getter]
    fn fingerprint(&self) -> String {
        self.inner.fingerprint().to_string()
    }

    #[getter]
    fn assembly(&self) -> PyAssemblyChoice {
        PyAssemblyChoice(self.inner.assembly.clone())
    }

    #[getter]
    fn model(&self) -> PyModelChoice {
        PyModelChoice(self.inner.model)
    }

    #[getter]
    fn altloc(&self) -> PyAltlocPolicy {
        PyAltlocPolicy(self.inner.altloc.clone())
    }

    #[getter]
    fn identifiers(&self) -> PyNamespace {
        self.inner.identifiers.into()
    }

    #[getter]
    fn missing_atoms(&self) -> PyMissingPolicy {
        self.inner.missing_atoms.into()
    }

    #[getter]
    fn hydrogens(&self) -> PyHydrogenPolicy {
        self.inner.hydrogens.into()
    }

    #[getter]
    fn atom_equivalence(&self) -> PyEquivalencePolicy {
        self.inner.atom_equivalence.into()
    }

    #[getter]
    fn symmetry(&self) -> PySymmetryPolicy {
        self.inner.symmetry.into()
    }

    #[getter]
    fn alignment(&self) -> PyAlignmentPolicy {
        self.inner.alignment.clone().into()
    }

    #[getter]
    fn precision(&self) -> PyPrecision {
        self.inner.precision.into()
    }

    #[getter]
    fn periodic(&self) -> PyPeriodicPolicy {
        self.inner.periodic.into()
    }

    #[getter]
    fn vdw_radii(&self) -> PyRadiusSet {
        match self.inner.vdw_radii {
            molframe::core::contract::RadiiSet::AmberUnited => PyRadiusSet::AmberUnited,
            molframe::core::contract::RadiiSet::Charmm => PyRadiusSet::Charmm,
            molframe::core::contract::RadiiSet::Alvarez => PyRadiusSet::Alvarez,
            _ => PyRadiusSet::Bondi,
        }
    }

    #[getter]
    fn contact_def(&self) -> PyContactDefinition {
        self.inner.contact_def.into()
    }

    #[getter]
    fn float_tolerance(&self) -> PyTolerance {
        self.inner.float_tolerance.into()
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

impl From<Namespace> for PyNamespace {
    fn from(value: Namespace) -> Self {
        match value {
            Namespace::Label => Self::Label,
            Namespace::Explicit => Self::Explicit,
            _ => Self::Auth,
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

impl From<MissingPolicy> for PyMissingPolicy {
    fn from(value: MissingPolicy) -> Self {
        match value {
            MissingPolicy::Ignore => Self::Ignore,
            MissingPolicy::Indeterminate => Self::Indeterminate,
            MissingPolicy::Fail => Self::Fail,
            _ => Self::Report,
        }
    }
}

impl From<ModelChoice> for PyModelChoice {
    fn from(value: ModelChoice) -> Self {
        Self(value)
    }
}

impl From<AltlocPolicy> for PyAltlocPolicy {
    fn from(value: AltlocPolicy) -> Self {
        Self(value)
    }
}

impl From<AssemblyChoice> for PyAssemblyChoice {
    fn from(value: AssemblyChoice) -> Self {
        Self(value)
    }
}

impl From<AnalysisPolicy> for PyAnalysisPolicy {
    fn from(inner: AnalysisPolicy) -> Self {
        Self { inner }
    }
}
