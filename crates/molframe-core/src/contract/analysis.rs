//! A result, and everything needed to judge it.
//!
//! An analysis returns a value *and* what it had to assume to produce one: how
//! much of the intended data was actually there, which decisions were made on
//! the caller's behalf, and where the inputs came from. The value alone is the
//! part every library gives you; the rest is what makes it citable.
//!
//! Declining to answer is a result in its own right. A library that cannot
//! return "no defensible number exists here" will always return a
//! defensible-looking wrong one instead. So an indeterminate analysis has no
//! value to misread: its [`Outcome`] holds the reason in the value's place.

use super::outcome::{Indeterminacy, Outcome, Quality};
use super::policy::{AnalysisPolicy, Fingerprint, PolicyField};
use super::provenance::Provenance;
use crate::diagnostic::Diagnostic;

/// How far an analysis got.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
#[non_exhaustive]
pub enum Status {
    /// Every intended atom was available and unambiguous.
    #[default]
    Complete,
    /// Computed, but some intended atoms were not available.
    Partial,
    /// Several defensible answers exist; one was chosen deterministically.
    Ambiguous,
    /// No defensible single answer exists under this policy.
    ///
    /// An analysis in this state has no value: see [`Outcome::Indeterminate`].
    Indeterminate,
}

impl Status {
    /// Returns true when the value may be used as an answer.
    #[must_use]
    pub const fn is_usable(self) -> bool {
        !matches!(self, Self::Indeterminate)
    }
}

/// How much of what an analysis meant to use it actually used.
///
/// "Intended" comes from what a component should contain, not from what the
/// file happened to have — so a residue whose side chain was never modelled
/// reduces coverage rather than silently shrinking the calculation.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Coverage {
    /// Atoms the analysis meant to use.
    pub intended: u32,
    /// Atoms it did use.
    pub used: u32,
    /// Atoms absent from the model.
    pub missing: u32,
    /// Atoms present but not uniquely resolvable.
    pub ambiguous: u32,
}

impl Coverage {
    /// Coverage for an analysis that used everything it meant to.
    #[must_use]
    pub const fn complete(atoms: u32) -> Self {
        Self {
            intended: atoms,
            used: atoms,
            missing: 0,
            ambiguous: 0,
        }
    }

    /// The fraction of intended atoms that were used.
    ///
    /// An analysis over two fifths of the intended atoms is not the same
    /// quantity as one over all of them, and this is what says so.
    ///
    /// # Panics
    /// The checked conversion accepts every unsigned 32-bit count.
    #[must_use]
    pub fn fraction(self) -> f64 {
        if self.intended == 0 {
            1.0
        } else {
            f64::from(self.used) / f64::from(self.intended)
        }
    }

    /// Returns true when nothing was missing or ambiguous.
    #[must_use]
    pub const fn is_complete(self) -> bool {
        self.missing == 0 && self.ambiguous == 0 && self.used == self.intended
    }
}

/// How much a decision changed the answer, where that is known.
///
/// An impact is a measurement or it is nothing. There is deliberately no
/// "low", "moderate" or "high": such a label would be a property an author
/// assigned, and presenting it beside measured values would lend it evidence it
/// does not have. The way to learn an impact is to vary the decision and look,
/// which is what an audit does and what [`MeasuredImpact`] records.
#[derive(Clone, Copy, PartialEq, Debug, Default)]
#[non_exhaustive]
pub enum Impact {
    /// Not measured. The honest default: nothing has been varied.
    #[default]
    Unmeasured,
    /// The decision cannot change this result, for a structural reason (the
    /// input has no alternate locations, say), not an estimated one.
    Inert,
    /// Measured by varying the decision across defensible alternatives.
    Measured(MeasuredImpact),
}

/// The measured effect of one decision on one analysis.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MeasuredImpact {
    /// The distance used, such as `"jaccard-distance"` or `"relative-error"`.
    pub metric: &'static str,
    /// Mean distance between results that differ in this decision alone, in the
    /// metric's own units.
    pub change: f64,
    /// How many pairs of results the measurement compared.
    pub comparisons: u32,
}

/// Where a decision came from.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[non_exhaustive]
pub enum AssumptionSource {
    /// The caller set it.
    Explicit,
    /// It came from the profile.
    ProfileDefault,
    /// It was inferred from the data.
    Inferred,
}

/// One decision made on the caller's behalf.
#[derive(Clone, PartialEq, Debug)]
pub struct Assumption {
    /// Which decision.
    pub field: PolicyField,
    /// What was decided.
    pub value: Box<str>,
    /// Where the decision came from.
    pub source: AssumptionSource,
    /// How much it changed the answer, if that has been measured.
    pub impact: Impact,
}

impl Assumption {
    /// Records a decision.
    #[must_use]
    pub fn new(
        field: PolicyField,
        value: impl Into<Box<str>>,
        source: AssumptionSource,
        impact: Impact,
    ) -> Self {
        Self {
            field,
            value: value.into(),
            source,
            impact,
        }
    }

    /// Returns true when the caller did not choose this and it measurably
    /// changed the answer by at least `threshold`, in the metric's units.
    ///
    /// These are the assumptions worth surfacing: a default the caller never saw
    /// that changed the answer is the failure this whole layer exists to expose.
    /// An unmeasured impact is never material by this test; it is unknown.
    #[must_use]
    pub fn is_silent_and_material(&self, threshold: f64) -> bool {
        self.source == AssumptionSource::ProfileDefault
            && matches!(self.impact, Impact::Measured(measured) if measured.change >= threshold)
    }
}

/// An outcome, and everything needed to judge it.
///
/// The answer is behind [`Analysis::outcome`], [`Analysis::value`] or
/// [`Analysis::into_result`], each of which says what happens when there is no
/// answer. There is no `Deref` to the value: a contract that lets the value be
/// read without asking whether it exists is a contract the indeterminate case
/// can slip through.
///
/// # Examples
///
/// ```
/// use molframe_core::contract::{Analysis, AnalysisPolicy, Coverage, Status};
///
/// let policy = AnalysisPolicy::default();
/// let radius = Analysis::complete(14.2_f64, Coverage::complete(1_960), &policy);
///
/// assert_eq!(radius.value(), Some(&14.2));
/// assert_eq!(radius.status(), Status::Complete);
/// assert_eq!(radius.coverage.fraction(), 1.0);
/// ```
#[derive(Clone, Debug)]
pub struct Analysis<T> {
    outcome: Outcome<T>,
    quality: Quality,
    /// How much of the intended data was used.
    pub coverage: Coverage,
    /// What the analysis found worth saying.
    pub warnings: Vec<Diagnostic>,
    /// What it decided on the caller's behalf.
    pub assumptions: Vec<Assumption>,
    /// Where the inputs came from and under what policy.
    pub provenance: Provenance,
}

impl<T> Analysis<T> {
    /// An answer of the given quality.
    #[must_use]
    pub fn determinate(
        value: T,
        quality: Quality,
        coverage: Coverage,
        policy: &AnalysisPolicy,
    ) -> Self {
        Self::from_parts(
            Outcome::Determinate(value),
            quality,
            coverage,
            Vec::new(),
            Vec::new(),
            Provenance::new(policy),
        )
    }

    /// A complete result.
    #[must_use]
    pub fn complete(value: T, coverage: Coverage, policy: &AnalysisPolicy) -> Self {
        Self::determinate(value, Quality::Complete, coverage, policy)
    }

    /// A result computed over less than it intended.
    #[must_use]
    pub fn partial(value: T, coverage: Coverage, policy: &AnalysisPolicy) -> Self {
        Self::determinate(value, Quality::Partial, coverage, policy)
    }

    /// A result that is one of several defensible ones.
    #[must_use]
    pub fn ambiguous(value: T, coverage: Coverage, policy: &AnalysisPolicy) -> Self {
        Self::determinate(value, Quality::Ambiguous, coverage, policy)
    }

    /// A refusal to answer, with the reason.
    ///
    /// There is no value to supply and none to read back.
    #[must_use]
    pub fn indeterminate(
        reason: Indeterminacy,
        coverage: Coverage,
        policy: &AnalysisPolicy,
    ) -> Self {
        Self::from_parts(
            Outcome::Indeterminate(reason),
            Quality::default(),
            coverage,
            Vec::new(),
            Vec::new(),
            Provenance::new(policy),
        )
    }

    /// Assembles an analysis from its parts.
    ///
    /// `quality` qualifies a determinate outcome and is ignored for an
    /// indeterminate one, which has no answer to qualify.
    #[must_use]
    pub const fn from_parts(
        outcome: Outcome<T>,
        quality: Quality,
        coverage: Coverage,
        warnings: Vec<Diagnostic>,
        assumptions: Vec<Assumption>,
        provenance: Provenance,
    ) -> Self {
        Self {
            outcome,
            quality,
            coverage,
            warnings,
            assumptions,
            provenance,
        }
    }

    /// Takes the analysis apart, for code that builds a new one from an old one.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        Outcome<T>,
        Quality,
        Coverage,
        Vec<Diagnostic>,
        Vec<Assumption>,
        Provenance,
    ) {
        (
            self.outcome,
            self.quality,
            self.coverage,
            self.warnings,
            self.assumptions,
            self.provenance,
        )
    }

    /// The answer or the reason there is none.
    #[must_use]
    pub const fn outcome(&self) -> &Outcome<T> {
        &self.outcome
    }

    /// Consumes the analysis, keeping only the outcome.
    #[must_use]
    pub fn into_outcome(self) -> Outcome<T> {
        self.outcome
    }

    /// The answer, if there is one.
    #[must_use]
    pub const fn value(&self) -> Option<&T> {
        self.outcome.value()
    }

    /// The reason there is no answer, if there is none.
    #[must_use]
    pub const fn indeterminacy(&self) -> Option<&Indeterminacy> {
        self.outcome.indeterminacy()
    }

    /// Consumes the analysis: the answer, or the reason it is missing.
    ///
    /// # Errors
    ///
    /// Returns the [`Indeterminacy`] when no defensible answer exists.
    pub fn into_result(self) -> Result<T, Indeterminacy> {
        self.outcome.into_result()
    }

    /// How far the analysis got: the quality of its answer, or `Indeterminate`.
    #[must_use]
    pub const fn status(&self) -> Status {
        match self.outcome {
            Outcome::Determinate(_) => match self.quality {
                Quality::Complete => Status::Complete,
                Quality::Partial => Status::Partial,
                Quality::Ambiguous => Status::Ambiguous,
            },
            Outcome::Indeterminate(_) => Status::Indeterminate,
        }
    }

    /// The quality of the answer, when there is one.
    #[must_use]
    pub const fn quality(&self) -> Option<Quality> {
        match self.outcome {
            Outcome::Determinate(_) => Some(self.quality),
            Outcome::Indeterminate(_) => None,
        }
    }

    /// Lowers the quality to `quality` if that is worse; no effect on a refusal.
    pub fn degrade(&mut self, quality: Quality) {
        self.quality = self.quality.worst(quality);
    }

    /// Attaches a finding.
    #[must_use]
    pub fn with_warning(mut self, warning: Diagnostic) -> Self {
        self.warnings.push(warning);
        self
    }

    /// Records a decision made on the caller's behalf.
    #[must_use]
    pub fn with_assumption(mut self, assumption: Assumption) -> Self {
        self.assumptions.push(assumption);
        self
    }

    /// The fingerprint of the policy this was computed under.
    #[must_use]
    pub fn policy_fingerprint(&self) -> Fingerprint {
        self.provenance.policy_fingerprint
    }

    /// The decisions the caller never made that measurably changed the answer by
    /// at least `threshold`.
    pub fn silent_assumptions(&self, threshold: f64) -> impl Iterator<Item = &Assumption> {
        self.assumptions
            .iter()
            .filter(move |assumption| assumption.is_silent_and_material(threshold))
    }

    /// Applies a function to the answer, keeping everything else.
    pub fn map<U>(self, transform: impl FnOnce(T) -> U) -> Analysis<U> {
        Analysis {
            outcome: self.outcome.map(transform),
            quality: self.quality,
            coverage: self.coverage,
            warnings: self.warnings,
            assumptions: self.assumptions,
            provenance: self.provenance,
        }
    }
}

#[cfg(test)]
#[path = "analysis_tests.rs"]
mod tests;
