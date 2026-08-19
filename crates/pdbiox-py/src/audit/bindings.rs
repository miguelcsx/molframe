//! Python projections for bounded policy audits.

#[path = "audit_plan.rs"]
mod audit_plan;
#[path = "audit_results.rs"]
mod audit_results;

pub(crate) use audit_plan::{PyAuditPlan, PyPolicySpace};

use crate::chemistry::PyRadiusSet;
use crate::query::{PyAltlocPolicy, PyAssemblyChoice, PyMissingPolicy, PyModelChoice, PyNamespace};
use pyo3::prelude::*;

macro_rules! policy_enum {
    ($python:literal, $name:ident, $rust:ty, {$($variant:ident),+ $(,)?}) => {
        #[pyclass(name = $python, frozen, eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) enum $name { $($variant),+ }
        impl From<$name> for $rust {
            fn from(value: $name) -> Self {
                match value { $( $name::$variant => <$rust>::$variant, )+ }
            }
        }
    };
}

policy_enum!("HydrogenPolicy", PyHydrogenPolicy, pdbiox::core::contract::HydrogenPolicy, {
    ExplicitOnly, IncludeInferred, Exclude
});
policy_enum!("EquivalencePolicy", PyEquivalencePolicy, pdbiox::core::contract::EquivalencePolicy, {
    None, Ccd, Explicit
});
policy_enum!("SymmetryPolicy", PySymmetryPolicy, pdbiox::core::contract::SymmetryPolicy, {
    None, Crystallographic, BiologicalAssembly
});
policy_enum!("Precision", PyPrecision, pdbiox::core::contract::Precision, { F32, F64 });
policy_enum!("PeriodicPolicy", PyPeriodicPolicy, pdbiox::core::contract::PeriodicPolicy, {
    None, Pbc, MinimumImage
});

impl From<pdbiox::core::contract::HydrogenPolicy> for PyHydrogenPolicy {
    fn from(value: pdbiox::core::contract::HydrogenPolicy) -> Self {
        match value {
            pdbiox::core::contract::HydrogenPolicy::IncludeInferred => Self::IncludeInferred,
            pdbiox::core::contract::HydrogenPolicy::Exclude => Self::Exclude,
            _ => Self::ExplicitOnly,
        }
    }
}

impl From<pdbiox::core::contract::EquivalencePolicy> for PyEquivalencePolicy {
    fn from(value: pdbiox::core::contract::EquivalencePolicy) -> Self {
        match value {
            pdbiox::core::contract::EquivalencePolicy::Explicit => Self::Explicit,
            _ => Self::Ccd,
        }
    }
}

impl From<pdbiox::core::contract::SymmetryPolicy> for PySymmetryPolicy {
    fn from(value: pdbiox::core::contract::SymmetryPolicy) -> Self {
        match value {
            pdbiox::core::contract::SymmetryPolicy::Crystallographic => Self::Crystallographic,
            pdbiox::core::contract::SymmetryPolicy::BiologicalAssembly => Self::BiologicalAssembly,
            _ => Self::None,
        }
    }
}

impl From<pdbiox::core::contract::Precision> for PyPrecision {
    fn from(value: pdbiox::core::contract::Precision) -> Self {
        match value {
            pdbiox::core::contract::Precision::F32 => Self::F32,
            _ => Self::F64,
        }
    }
}

impl From<pdbiox::core::contract::PeriodicPolicy> for PyPeriodicPolicy {
    fn from(value: pdbiox::core::contract::PeriodicPolicy) -> Self {
        match value {
            pdbiox::core::contract::PeriodicPolicy::Pbc => Self::Pbc,
            pdbiox::core::contract::PeriodicPolicy::MinimumImage => Self::MinimumImage,
            _ => Self::None,
        }
    }
}

impl From<pdbiox::core::contract::AlignmentPolicy> for PyAlignmentPolicy {
    fn from(value: pdbiox::core::contract::AlignmentPolicy) -> Self {
        match value {
            pdbiox::core::contract::AlignmentPolicy::Explicit(selection) => {
                Self::Explicit(selection.into())
            }
            pdbiox::core::contract::AlignmentPolicy::Global => Self::Global(),
            pdbiox::core::contract::AlignmentPolicy::Local => Self::Local(),
            _ => Self::None(),
        }
    }
}

impl From<pdbiox::core::contract::ContactDefinition> for PyContactDefinition {
    fn from(value: pdbiox::core::contract::ContactDefinition) -> Self {
        match value {
            pdbiox::core::contract::ContactDefinition::DistanceCutoff { tolerance } => {
                Self::DistanceCutoff { tolerance }
            }
            pdbiox::core::contract::ContactDefinition::SurfaceBased { probe } => {
                Self::SurfaceBased { probe }
            }
            _ => Self::DistanceCutoff { tolerance: 0.5 },
        }
    }
}

#[pymethods]
impl PyEquivalencePolicy {
    #[staticmethod]
    fn none() -> Self {
        Self::None
    }
}

#[pymethods]
impl PySymmetryPolicy {
    #[staticmethod]
    fn none() -> Self {
        Self::None
    }
}

#[pymethods]
impl PyPeriodicPolicy {
    #[staticmethod]
    fn none() -> Self {
        Self::None
    }
}

#[pyclass(name = "AlignmentPolicy", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) enum PyAlignmentPolicy {
    None(),
    Explicit(String),
    Global(),
    Local(),
}

impl From<PyAlignmentPolicy> for pdbiox::core::contract::AlignmentPolicy {
    fn from(value: PyAlignmentPolicy) -> Self {
        match value {
            PyAlignmentPolicy::None() => Self::None,
            PyAlignmentPolicy::Explicit(value) => Self::Explicit(value.into()),
            PyAlignmentPolicy::Global() => Self::Global,
            PyAlignmentPolicy::Local() => Self::Local,
        }
    }
}

#[pymethods]
impl PyAlignmentPolicy {
    #[staticmethod]
    fn none() -> Self {
        Self::None()
    }
    #[staticmethod]
    fn explicit(selection: String) -> Self {
        Self::Explicit(selection)
    }
    #[staticmethod]
    #[pyo3(name = "global_")]
    fn global_() -> Self {
        Self::Global()
    }
    #[staticmethod]
    fn local() -> Self {
        Self::Local()
    }
}

#[pyclass(name = "ContactDefinition", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) enum PyContactDefinition {
    DistanceCutoff { tolerance: f32 },
    SurfaceBased { probe: f32 },
}

impl From<PyContactDefinition> for pdbiox::core::contract::ContactDefinition {
    fn from(value: PyContactDefinition) -> Self {
        match value {
            PyContactDefinition::DistanceCutoff { tolerance } => Self::DistanceCutoff { tolerance },
            PyContactDefinition::SurfaceBased { probe } => Self::SurfaceBased { probe },
        }
    }
}

#[pymethods]
impl PyContactDefinition {
    #[staticmethod]
    fn distance_cutoff(tolerance: f32) -> Self {
        Self::DistanceCutoff { tolerance }
    }
    #[staticmethod]
    fn surface_based(probe: f32) -> Self {
        Self::SurfaceBased { probe }
    }
}

#[pyclass(name = "Tolerance", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyTolerance {
    #[pyo3(get)]
    relative: f64,
    #[pyo3(get)]
    absolute: f64,
}

impl From<PyTolerance> for pdbiox::core::contract::Tolerance {
    fn from(value: PyTolerance) -> Self {
        Self {
            relative: value.relative,
            absolute: value.absolute,
        }
    }
}

impl From<pdbiox::core::contract::Tolerance> for PyTolerance {
    fn from(value: pdbiox::core::contract::Tolerance) -> Self {
        Self {
            relative: value.relative,
            absolute: value.absolute,
        }
    }
}

#[pymethods]
impl PyTolerance {
    #[new]
    fn new(relative: f64, absolute: f64) -> Self {
        Self { relative, absolute }
    }
    #[staticmethod]
    fn standard() -> Self {
        <Self as From<pdbiox::core::contract::Tolerance>>::from(
            pdbiox::core::contract::Tolerance::default(),
        )
    }
}

#[pyclass(name = "PolicyField", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPolicyField {
    Assembly,
    Model,
    Altloc,
    Identifiers,
    MissingAtoms,
    Hydrogens,
    AtomEquivalence,
    Symmetry,
    Alignment,
    Precision,
    Periodic,
    VdwRadii,
    ContactDef,
    FloatTolerance,
}

impl From<PyPolicyField> for pdbiox::core::contract::PolicyField {
    fn from(value: PyPolicyField) -> Self {
        match value {
            PyPolicyField::Assembly => Self::Assembly,
            PyPolicyField::Model => Self::Model,
            PyPolicyField::Altloc => Self::Altloc,
            PyPolicyField::Identifiers => Self::Identifiers,
            PyPolicyField::MissingAtoms => Self::MissingAtoms,
            PyPolicyField::Hydrogens => Self::Hydrogens,
            PyPolicyField::AtomEquivalence => Self::AtomEquivalence,
            PyPolicyField::Symmetry => Self::Symmetry,
            PyPolicyField::Alignment => Self::Alignment,
            PyPolicyField::Precision => Self::Precision,
            PyPolicyField::Periodic => Self::Periodic,
            PyPolicyField::VdwRadii => Self::VdwRadii,
            PyPolicyField::ContactDef => Self::ContactDef,
            PyPolicyField::FloatTolerance => Self::FloatTolerance,
        }
    }
}

impl From<pdbiox::core::contract::PolicyField> for PyPolicyField {
    fn from(value: pdbiox::core::contract::PolicyField) -> Self {
        match value {
            pdbiox::core::contract::PolicyField::Assembly => Self::Assembly,
            pdbiox::core::contract::PolicyField::Model => Self::Model,
            pdbiox::core::contract::PolicyField::Altloc => Self::Altloc,
            pdbiox::core::contract::PolicyField::Identifiers => Self::Identifiers,
            pdbiox::core::contract::PolicyField::MissingAtoms => Self::MissingAtoms,
            pdbiox::core::contract::PolicyField::Hydrogens => Self::Hydrogens,
            pdbiox::core::contract::PolicyField::AtomEquivalence => Self::AtomEquivalence,
            pdbiox::core::contract::PolicyField::Symmetry => Self::Symmetry,
            pdbiox::core::contract::PolicyField::Alignment => Self::Alignment,
            pdbiox::core::contract::PolicyField::Precision => Self::Precision,
            pdbiox::core::contract::PolicyField::Periodic => Self::Periodic,
            pdbiox::core::contract::PolicyField::VdwRadii => Self::VdwRadii,
            pdbiox::core::contract::PolicyField::ContactDef => Self::ContactDef,
            _ => Self::FloatTolerance,
        }
    }
}

#[pyclass(name = "PolicyValue", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolicyValue {
    #[pyo3(get)]
    field: PyPolicyField,
    #[pyo3(get)]
    value: String,
}

#[pymethods]
impl PyPolicyValue {
    #[staticmethod]
    fn named(field: PyPolicyField, value: String) -> Self {
        Self { field, value }
    }
}

#[pyclass(name = "PolicyDimension", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolicyDimension(pub(crate) pdbiox::PolicyDimension);

#[pymethods]
impl PyPolicyDimension {
    #[staticmethod]
    fn assembly(values: Vec<PyAssemblyChoice>) -> Self {
        Self(pdbiox::PolicyDimension::assembly(
            values.into_iter().map(|value| value.0),
        ))
    }
    #[staticmethod]
    fn model(values: Vec<PyModelChoice>) -> Self {
        Self(pdbiox::PolicyDimension::model(
            values.into_iter().map(|value| value.0),
        ))
    }
    #[staticmethod]
    fn altloc(values: Vec<PyAltlocPolicy>) -> Self {
        Self(pdbiox::PolicyDimension::altloc(
            values.into_iter().map(|value| value.0),
        ))
    }
    #[staticmethod]
    fn identifiers(values: Vec<PyNamespace>) -> Self {
        Self(pdbiox::PolicyDimension::identifiers(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn missing_atoms(values: Vec<PyMissingPolicy>) -> Self {
        Self(pdbiox::PolicyDimension::missing_atoms(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn hydrogens(values: Vec<PyHydrogenPolicy>) -> Self {
        Self(pdbiox::PolicyDimension::hydrogens(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn atom_equivalence(values: Vec<PyEquivalencePolicy>) -> Self {
        Self(pdbiox::PolicyDimension::atom_equivalence(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn symmetry(values: Vec<PySymmetryPolicy>) -> Self {
        Self(pdbiox::PolicyDimension::symmetry(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn alignment(values: Vec<PyAlignmentPolicy>) -> Self {
        Self(pdbiox::PolicyDimension::alignment(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn precision(values: Vec<PyPrecision>) -> Self {
        Self(pdbiox::PolicyDimension::precision(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn periodic(values: Vec<PyPeriodicPolicy>) -> Self {
        Self(pdbiox::PolicyDimension::periodic(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn vdw_radii(values: Vec<PyRadiusSet>) -> Self {
        Self(pdbiox::PolicyDimension::vdw_radii(values.into_iter().map(
            |value| match value {
                PyRadiusSet::Bondi => pdbiox::core::contract::RadiiSet::Bondi,
                PyRadiusSet::AmberUnited => pdbiox::core::contract::RadiiSet::AmberUnited,
                PyRadiusSet::Charmm => pdbiox::core::contract::RadiiSet::Charmm,
                PyRadiusSet::Alvarez => pdbiox::core::contract::RadiiSet::Alvarez,
            },
        )))
    }
    #[staticmethod]
    fn contact_def(values: Vec<PyContactDefinition>) -> Self {
        Self(pdbiox::PolicyDimension::contact_def(
            values.into_iter().map(Into::into),
        ))
    }
    #[staticmethod]
    fn float_tolerance(values: Vec<PyTolerance>) -> Self {
        Self(pdbiox::PolicyDimension::float_tolerance(
            values.into_iter().map(Into::into),
        ))
    }
    #[getter]
    fn field(&self) -> PyPolicyField {
        self.0.field().into()
    }
    fn __len__(&self) -> usize {
        self.0.len()
    }
    #[getter]
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

#[pyclass(name = "PlanError", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPlanError {
    #[pyo3(get)]
    message: String,
    #[pyo3(get)]
    field: Option<PyPolicyField>,
    #[pyo3(get)]
    cost: Option<usize>,
    #[pyo3(get)]
    limit: Option<usize>,
}

impl From<pdbiox::PlanError> for PyPlanError {
    fn from(value: pdbiox::PlanError) -> Self {
        let (field, cost, limit) = match value {
            pdbiox::PlanError::EmptyDimension(field)
            | pdbiox::PlanError::DuplicateDimension(field) => (Some(field.into()), None, None),
            pdbiox::PlanError::LimitExceeded { cost, limit } => (None, Some(cost), Some(limit)),
            _ => (None, None, None),
        };
        Self {
            message: value.to_string(),
            field,
            cost,
            limit,
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyHydrogenPolicy>()?;
    module.add_class::<PyEquivalencePolicy>()?;
    module.add_class::<PySymmetryPolicy>()?;
    module.add_class::<PyPrecision>()?;
    module.add_class::<PyPeriodicPolicy>()?;
    module.add_class::<PyAlignmentPolicy>()?;
    module.add_class::<PyContactDefinition>()?;
    module.add_class::<PyTolerance>()?;
    module.add_class::<PyPolicyField>()?;
    module.add_class::<PyPolicyValue>()?;
    module.add_class::<PyPolicyDimension>()?;
    module.add_class::<PyPlanError>()?;
    module.add_class::<PyPolicySpace>()?;
    module.add_class::<PyAuditPlan>()?;
    audit_results::register(module)?;
    Ok(())
}

#[cfg(test)]
#[path = "audit_tests.rs"]
mod tests;
