use crate::PolicyValue;
use pdbiox_core::contract::{
    AlignmentPolicy, AltlocPolicy, AnalysisPolicy, AssemblyChoice, ContactDefinition,
    EquivalencePolicy, HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy,
    PolicyField, Precision, RadiiSet, SymmetryPolicy, Tolerance,
};
use std::fmt;

/// One named policy field and the defensible values it may take.
#[derive(Clone, Debug, PartialEq)]
pub struct PolicyDimension {
    field: PolicyField,
    values: Vec<PolicyValue>,
}

macro_rules! dimension_constructor {
    ($name:ident, $field:ident, $variant:ident, $ty:ty) => {
        #[doc = concat!("Varies `", stringify!($field), "` over the supplied values.")]
        pub fn $name(values: impl IntoIterator<Item = $ty>) -> Self {
            Self {
                field: PolicyField::$field,
                values: values.into_iter().map(PolicyValue::$variant).collect(),
            }
        }
    };
}

impl PolicyDimension {
    dimension_constructor!(assembly, Assembly, Assembly, AssemblyChoice);
    dimension_constructor!(model, Model, Model, ModelChoice);
    dimension_constructor!(altloc, Altloc, Altloc, AltlocPolicy);
    dimension_constructor!(identifiers, Identifiers, Identifiers, Namespace);
    dimension_constructor!(missing_atoms, MissingAtoms, MissingAtoms, MissingPolicy);
    dimension_constructor!(hydrogens, Hydrogens, Hydrogens, HydrogenPolicy);
    dimension_constructor!(
        atom_equivalence,
        AtomEquivalence,
        AtomEquivalence,
        EquivalencePolicy
    );
    dimension_constructor!(symmetry, Symmetry, Symmetry, SymmetryPolicy);
    dimension_constructor!(alignment, Alignment, Alignment, AlignmentPolicy);
    dimension_constructor!(precision, Precision, Precision, Precision);
    dimension_constructor!(periodic, Periodic, Periodic, PeriodicPolicy);
    dimension_constructor!(vdw_radii, VdwRadii, VdwRadii, RadiiSet);
    dimension_constructor!(contact_def, ContactDef, ContactDef, ContactDefinition);
    dimension_constructor!(float_tolerance, FloatTolerance, FloatTolerance, Tolerance);

    /// The policy field varied by this dimension.
    #[must_use]
    pub const fn field(&self) -> PolicyField {
        self.field
    }

    /// The number of alternatives in this dimension.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether this dimension has no alternatives.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

/// A bounded Cartesian policy space that has not yet been expanded.
#[derive(Clone, Debug, PartialEq)]
pub struct PolicySpace {
    baseline: AnalysisPolicy,
    dimensions: Vec<PolicyDimension>,
    max_runs: usize,
}

impl PolicySpace {
    /// Starts a policy space from a baseline policy.
    #[must_use]
    pub const fn new(baseline: AnalysisPolicy) -> Self {
        Self {
            baseline,
            dimensions: Vec::new(),
            max_runs: 4_096,
        }
    }

    /// Adds a policy dimension.
    #[must_use]
    pub fn vary(mut self, dimension: PolicyDimension) -> Self {
        self.dimensions.push(dimension);
        self
    }

    /// Sets the hard upper bound on executed policy points.
    #[must_use]
    pub const fn with_max_runs(mut self, max_runs: usize) -> Self {
        self.max_runs = max_runs;
        self
    }

    /// Computes the exact number of policy points without expanding them.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty or duplicate dimension, or if the
    /// Cartesian product cannot be represented by `usize`.
    pub fn cost(&self) -> Result<usize, PlanError> {
        validate_dimensions(&self.dimensions)?;
        let mut cost = 1usize;
        for dimension in &self.dimensions {
            cost = cost
                .checked_mul(dimension.len())
                .ok_or(PlanError::CostOverflow)?;
        }
        Ok(cost)
    }

    /// Validates and deterministically expands the Cartesian policy space.
    ///
    /// # Errors
    ///
    /// Returns the validation errors from [`Self::cost`] or
    /// [`PlanError::LimitExceeded`] when the configured bound is too small.
    pub fn plan(self) -> Result<AuditPlan, PlanError> {
        let cost = self.cost()?;
        if cost > self.max_runs {
            return Err(PlanError::LimitExceeded {
                cost,
                limit: self.max_runs,
            });
        }
        let mut policies = vec![self.baseline];
        let mut coordinates = vec![Vec::new()];
        for dimension in &self.dimensions {
            let next_policy_count = policies
                .len()
                .checked_mul(dimension.len())
                .ok_or(PlanError::CostOverflow)?;
            let next_coordinate_count = coordinates
                .len()
                .checked_mul(dimension.len())
                .ok_or(PlanError::CostOverflow)?;
            let mut next_policies = Vec::with_capacity(next_policy_count);
            let mut next_coordinates = Vec::with_capacity(next_coordinate_count);
            for (policy, coordinate) in policies.iter().zip(&coordinates) {
                for (choice, value) in dimension.values.iter().enumerate() {
                    let mut varied = policy.clone();
                    apply(&mut varied, value);
                    let mut point = coordinate.clone();
                    point.push(choice);
                    next_policies.push(varied);
                    next_coordinates.push(point);
                }
            }
            policies = next_policies;
            coordinates = next_coordinates;
        }
        Ok(AuditPlan {
            fields: self.dimensions.iter().map(PolicyDimension::field).collect(),
            policies,
            coordinates,
        })
    }
}

fn validate_dimensions(dimensions: &[PolicyDimension]) -> Result<(), PlanError> {
    let mut fields = Vec::with_capacity(dimensions.len());
    for dimension in dimensions {
        if dimension.is_empty() {
            return Err(PlanError::EmptyDimension(dimension.field()));
        }
        if fields.contains(&dimension.field()) {
            return Err(PlanError::DuplicateDimension(dimension.field()));
        }
        fields.push(dimension.field());
    }
    Ok(())
}

fn apply(policy: &mut AnalysisPolicy, value: &PolicyValue) {
    match value {
        PolicyValue::Assembly(value) => policy.assembly = value.clone(),
        PolicyValue::Model(value) => policy.model = *value,
        PolicyValue::Altloc(value) => policy.altloc = value.clone(),
        PolicyValue::Identifiers(value) => policy.identifiers = *value,
        PolicyValue::MissingAtoms(value) => policy.missing_atoms = *value,
        PolicyValue::Hydrogens(value) => policy.hydrogens = *value,
        PolicyValue::AtomEquivalence(value) => policy.atom_equivalence = *value,
        PolicyValue::Symmetry(value) => policy.symmetry = *value,
        PolicyValue::Alignment(value) => policy.alignment = value.clone(),
        PolicyValue::Precision(value) => policy.precision = *value,
        PolicyValue::Periodic(value) => policy.periodic = *value,
        PolicyValue::VdwRadii(value) => policy.vdw_radii = *value,
        PolicyValue::ContactDef(value) => policy.contact_def = *value,
        PolicyValue::FloatTolerance(value) => policy.float_tolerance = *value,
    }
}

/// A validated, expanded audit plan.
#[derive(Clone, Debug, PartialEq)]
pub struct AuditPlan {
    pub(crate) fields: Vec<PolicyField>,
    pub(crate) policies: Vec<AnalysisPolicy>,
    pub(crate) coordinates: Vec<Vec<usize>>,
}

impl AuditPlan {
    /// Exact number of analysis runs.
    #[must_use]
    pub fn cost(&self) -> usize {
        self.policies.len()
    }

    /// Fields varied by the plan, in declaration order.
    #[must_use]
    pub fn fields(&self) -> &[PolicyField] {
        &self.fields
    }

    /// Expanded policies, in deterministic Cartesian order.
    #[must_use]
    pub fn policies(&self) -> &[AnalysisPolicy] {
        &self.policies
    }

    /// Cartesian coordinates identifying the selected alternative in each field.
    #[must_use]
    pub fn coordinates(&self) -> &[Vec<usize>] {
        &self.coordinates
    }
}

/// Why a policy space could not be planned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum PlanError {
    /// A named dimension has no alternatives.
    EmptyDimension(PolicyField),
    /// A field was declared more than once.
    DuplicateDimension(PolicyField),
    /// The Cartesian product overflowed the platform's index size.
    CostOverflow,
    /// The requested space exceeds the caller's bound.
    LimitExceeded {
        /// Exact requested number of runs.
        cost: usize,
        /// Configured upper bound.
        limit: usize,
    },
}

impl fmt::Display for PlanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDimension(field) => {
                write!(formatter, "{} has no alternatives", field.name())
            }
            Self::DuplicateDimension(field) => {
                write!(formatter, "{} is varied twice", field.name())
            }
            Self::CostOverflow => formatter.write_str("policy-space cost overflowed usize"),
            Self::LimitExceeded { cost, limit } => {
                write!(
                    formatter,
                    "policy space needs {cost} runs, exceeding limit {limit}"
                )
            }
        }
    }
}

impl std::error::Error for PlanError {}

#[cfg(test)]
#[path = "plan_tests.rs"]
mod tests;
