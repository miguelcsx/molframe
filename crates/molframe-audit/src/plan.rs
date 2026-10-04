use crate::{PlanError, PolicyValue, UncertaintyClass};
use molframe_core::contract::{
    AlignmentPolicy, AltlocPolicy, AnalysisPolicy, AssemblyChoice, ContactDefinition,
    EquivalencePolicy, HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy,
    PolicyField, PolicyParseError, Precision, RadiiSet, SymmetryPolicy, Tolerance,
};

/// One named policy field and the defensible values it may take.
#[derive(Clone, Debug, PartialEq)]
pub struct PolicyDimension {
    field: PolicyField,
    values: Vec<PolicyValue>,
    rationale: Box<str>,
    evidence: Box<str>,
}

macro_rules! dimension_constructor {
    ($name:ident, $field:ident, $variant:ident, $ty:ty) => {
        #[doc = concat!("Varies `", stringify!($field), "` over the supplied values.")]
        pub fn $name(values: impl IntoIterator<Item = $ty>) -> Self {
            Self {
                field: PolicyField::$field,
                values: values.into_iter().map(PolicyValue::$variant).collect(),
                rationale: Box::default(),
                evidence: Box::default(),
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

    /// Varies `field` over alternatives written in the vocabulary a written policy uses.
    ///
    /// Each word is read by the same parser as the policy itself, so `"first"`,
    /// `"biological:1"` and `"crystal:8.0"` mean here what they mean there.
    /// `float_tolerance` alternatives are written `"relative,absolute"`.
    ///
    /// # Errors
    ///
    /// Returns the policy parse error for the first word the field's vocabulary
    /// does not contain.
    pub fn named(field: PolicyField, words: &[&str]) -> Result<Self, PolicyParseError> {
        fn all<T: std::str::FromStr<Err = PolicyParseError>>(
            words: &[&str],
            wrap: fn(T) -> PolicyValue,
        ) -> Result<Vec<PolicyValue>, PolicyParseError> {
            words.iter().map(|word| word.parse().map(wrap)).collect()
        }
        let values = match field {
            PolicyField::Assembly => all(words, PolicyValue::Assembly)?,
            PolicyField::Model => all(words, PolicyValue::Model)?,
            PolicyField::Altloc => all(words, PolicyValue::Altloc)?,
            PolicyField::Identifiers => all(words, PolicyValue::Identifiers)?,
            PolicyField::MissingAtoms => all(words, PolicyValue::MissingAtoms)?,
            PolicyField::Hydrogens => all(words, PolicyValue::Hydrogens)?,
            PolicyField::AtomEquivalence => all(words, PolicyValue::AtomEquivalence)?,
            PolicyField::Symmetry => all(words, PolicyValue::Symmetry)?,
            PolicyField::Alignment => all(words, PolicyValue::Alignment)?,
            PolicyField::Precision => all(words, PolicyValue::Precision)?,
            PolicyField::Periodic => all(words, PolicyValue::Periodic)?,
            PolicyField::VdwRadii => all(words, PolicyValue::VdwRadii)?,
            PolicyField::ContactDef => all(words, PolicyValue::ContactDef)?,
            PolicyField::FloatTolerance => words
                .iter()
                .map(|word| tolerance(word).map(PolicyValue::FloatTolerance))
                .collect::<Result<_, _>>()?,
            _ => {
                return Err(PolicyParseError::new(
                    "policy field",
                    field.name(),
                    "a field this version of the audit knows",
                ));
            }
        };
        Ok(Self {
            field,
            values,
            rationale: Box::default(),
            evidence: Box::default(),
        })
    }

    /// Records why these alternatives are the defensible ones, and what supports that.
    ///
    /// A sweep over alternatives nobody can defend measures nothing about the
    /// analysis. The rationale and the evidence travel with the plan, so a reader of
    /// a result can see what was claimed to be reasonable and on what ground.
    #[must_use]
    pub fn justified(
        mut self,
        rationale: impl Into<Box<str>>,
        evidence: impl Into<Box<str>>,
    ) -> Self {
        self.rationale = rationale.into();
        self.evidence = evidence.into();
        self
    }

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
    exclusions: Vec<(PolicyValue, PolicyValue)>,
    max_runs: usize,
}

impl PolicySpace {
    /// Starts a policy space from a baseline policy.
    #[must_use]
    pub const fn new(baseline: AnalysisPolicy) -> Self {
        Self {
            baseline,
            dimensions: Vec::new(),
            exclusions: Vec::new(),
            max_runs: 4_096,
        }
    }

    /// Adds a policy dimension.
    #[must_use]
    pub fn vary(mut self, dimension: PolicyDimension) -> Self {
        self.dimensions.push(dimension);
        self
    }

    /// Declares that a universe holding `when` cannot also hold `then_not`.
    ///
    /// Some pairs of decisions are not both defensible, whatever the policy
    /// consistency check says: a rule that discards altlocs makes an occupancy sum
    /// meaningless. Such a pair is removed from a constrained plan and refuses a strict one.
    #[must_use]
    pub fn forbid(mut self, when: PolicyValue, then_not: PolicyValue) -> Self {
        self.exclusions.push((when, then_not));
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
    /// The plan is the whole product, so every main effect and interaction is
    /// defined. A combination that cannot be run, or that the space forbids, is
    /// refused rather than dropped: dropping it would leave a plan that is not the
    /// product and that the additive split cannot describe. Use
    /// [`Self::plan_constrained`] to drop such combinations deliberately.
    ///
    /// # Errors
    ///
    /// Returns the validation errors from [`Self::cost`], [`PlanError::LimitExceeded`]
    /// when the configured bound is too small, and [`PlanError::Conflict`] for a
    /// combination that contradicts itself or is forbidden.
    pub fn plan(self) -> Result<AuditPlan, PlanError> {
        let (policies, coordinates) = self.expand()?;
        if let Some(run) = policies.iter().position(|policy| !self.allowed(policy)) {
            return Err(PlanError::Conflict { run });
        }
        Ok(self.finish(policies, coordinates, 0))
    }

    /// Expands the space and keeps only the combinations that can be run and are not forbidden.
    ///
    /// The plan is then not the whole product: the additive main-effect and interaction
    /// shares do not apply, and the Shapley shares, which do not need a product, carry the
    /// attribution. The number of combinations dropped is recorded.
    ///
    /// # Errors
    ///
    /// Returns the validation errors from [`Self::cost`], [`PlanError::LimitExceeded`] when
    /// the bound is too small, and [`PlanError::NoUniverse`] when nothing survives.
    pub fn plan_constrained(self) -> Result<AuditPlan, PlanError> {
        let (policies, coordinates) = self.expand()?;
        let total = policies.len();
        let (policies, coordinates): (Vec<_>, Vec<_>) = policies
            .into_iter()
            .zip(coordinates)
            .filter(|(policy, _)| self.allowed(policy))
            .unzip();
        if policies.is_empty() {
            return Err(PlanError::NoUniverse);
        }
        let skipped = total - policies.len();
        Ok(self.finish(policies, coordinates, skipped))
    }

    fn allowed(&self, policy: &AnalysisPolicy) -> bool {
        policy.check_consistency().is_ok()
            && !self
                .exclusions
                .iter()
                .any(|(when, then_not)| holds(policy, when) && holds(policy, then_not))
    }

    fn finish(
        self,
        policies: Vec<AnalysisPolicy>,
        coordinates: Vec<Vec<usize>>,
        skipped: usize,
    ) -> AuditPlan {
        let levels_total: usize = self.dimensions.iter().map(PolicyDimension::len).product();
        AuditPlan {
            fields: self.dimensions.iter().map(PolicyDimension::field).collect(),
            balanced: policies.len() == levels_total,
            skipped,
            decisions: self
                .dimensions
                .iter()
                .map(|dimension| Decision {
                    field: dimension.field,
                    class: UncertaintyClass::of(dimension.field),
                    rationale: dimension.rationale.clone(),
                    evidence: dimension.evidence.clone(),
                })
                .collect(),
            policies,
            coordinates,
        }
    }

    fn expand(&self) -> Result<(Vec<AnalysisPolicy>, Vec<Vec<usize>>), PlanError> {
        let cost = self.cost()?;
        if cost > self.max_runs {
            return Err(PlanError::LimitExceeded {
                cost,
                limit: self.max_runs,
            });
        }
        let mut policies = vec![self.baseline.clone()];
        let mut coordinates = vec![Vec::new()];
        for dimension in &self.dimensions {
            let mut next_policies = Vec::with_capacity(policies.len() * dimension.len());
            let mut next_coordinates = Vec::with_capacity(coordinates.len() * dimension.len());
            for (policy, coordinate) in policies.iter().zip(&coordinates) {
                for (choice, value) in dimension.values.iter().enumerate() {
                    let mut varied = policy.clone();
                    apply(&mut varied, value);
                    let mut point: Vec<usize> = coordinate.clone();
                    point.push(choice);
                    next_policies.push(varied);
                    next_coordinates.push(point);
                }
            }
            policies = next_policies;
            coordinates = next_coordinates;
        }
        Ok((policies, coordinates))
    }
}

fn tolerance(word: &str) -> Result<Tolerance, PolicyParseError> {
    let refuse = || PolicyParseError::new("float_tolerance", word, "<relative>,<absolute>");
    let (relative, absolute) = word.split_once(',').ok_or_else(refuse)?;
    let parse = |text: &str| {
        text.trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v >= 0.0)
    };
    match (parse(relative), parse(absolute)) {
        (Some(relative), Some(absolute)) => Ok(Tolerance { relative, absolute }),
        _ => Err(refuse()),
    }
}

/// Whether the policy holds `value` for the field it names.
fn holds(policy: &AnalysisPolicy, value: &PolicyValue) -> bool {
    let mut probe = policy.clone();
    apply(&mut probe, value);
    probe == *policy
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
    pub(crate) balanced: bool,
    pub(crate) skipped: usize,
    pub(crate) decisions: Vec<Decision>,
}

/// A varied decision as the plan records it: its class, and why it is defensible to vary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decision {
    /// The decision.
    pub field: PolicyField,
    /// The kind of uncertainty it carries.
    pub class: UncertaintyClass,
    /// Why the alternatives are the defensible ones.
    pub rationale: Box<str>,
    /// What supports that.
    pub evidence: Box<str>,
}

impl AuditPlan {
    /// Exact number of analysis runs.
    #[must_use]
    pub fn cost(&self) -> usize {
        self.policies.len()
    }

    /// Whether the plan is the whole product of the alternatives.
    ///
    /// Only then do the main-effect and interaction shares of a decomposition add up.
    #[must_use]
    pub const fn balanced(&self) -> bool {
        self.balanced
    }

    /// How many combinations a constrained plan dropped.
    #[must_use]
    pub const fn skipped(&self) -> usize {
        self.skipped
    }

    /// The varied decisions with their class and justification.
    #[must_use]
    pub fn decisions(&self) -> &[Decision] {
        &self.decisions
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

    /// Refuses a plan that varies a field the analysis does not apply.
    ///
    /// `read` is what the analysis recorded as applied. A decision the analysis
    /// never sees cannot move its answer, so a sweep over it would come back
    /// perfectly stable and mean nothing.
    ///
    /// # Errors
    ///
    /// Returns [`PlanError::NotRead`] for the first varied field not in `read`.
    pub fn require_read(&self, read: &[PolicyField]) -> Result<(), PlanError> {
        match self.fields.iter().find(|field| !read.contains(field)) {
            Some(field) => Err(PlanError::NotRead(*field)),
            None => Ok(()),
        }
    }

    /// Cartesian coordinates identifying the selected alternative in each field.
    #[must_use]
    pub fn coordinates(&self) -> &[Vec<usize>] {
        &self.coordinates
    }
}

#[cfg(test)]
#[path = "plan_tests.rs"]
mod tests;
