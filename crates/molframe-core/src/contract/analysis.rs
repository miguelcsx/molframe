//! A result, and everything needed to judge it.
//!
//! An analysis returns a value *and* what it had to assume to produce one: how
//! much of the intended data was actually there, which decisions were made on
//! the caller's behalf, and where the inputs came from. The value alone is the
//! part every library gives you; the rest is what makes it citable.
//!
//! Declining to answer is a result in its own right. A library that cannot
//! return "no defensible number exists here" will always return a
//! defensible-looking wrong one instead.

use super::policy::{AnalysisPolicy, Fingerprint, PolicyField};
use super::provenance::Provenance;
use crate::diagnostic::Diagnostic;
use std::ops::Deref;

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
    /// The value carried alongside this status is not meaningful.
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

/// How much a decision could have changed the answer.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
#[non_exhaustive]
pub enum ImpactEstimate {
    /// The decision could not have changed anything here.
    None,
    /// Unlikely to matter.
    Low,
    /// Could matter.
    Moderate,
    /// Likely to matter.
    High,
    /// Not estimated.
    #[default]
    Unknown,
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
    /// How much it could have mattered.
    pub impact: ImpactEstimate,
}

impl Assumption {
    /// Records a decision.
    #[must_use]
    pub fn new(
        field: PolicyField,
        value: impl Into<Box<str>>,
        source: AssumptionSource,
        impact: ImpactEstimate,
    ) -> Self {
        Self {
            field,
            value: value.into(),
            source,
            impact,
        }
    }

    /// Returns true when the caller did not choose this and it could matter.
    ///
    /// These are the assumptions worth surfacing: a default the caller never saw
    /// that changed the answer is the failure this whole layer exists to expose.
    #[must_use]
    pub fn is_silent_and_material(&self) -> bool {
        self.source == AssumptionSource::ProfileDefault
            && self.impact >= ImpactEstimate::Moderate
            && self.impact != ImpactEstimate::Unknown
    }
}

/// A value, and everything needed to judge it.
///
/// Dereferences to the value, so reading the answer costs one character and
/// reading the rest is there when it matters. Ergonomics are load-bearing here:
/// a contract that makes ordinary work painful gets bypassed, and a bypassed
/// contract records nothing.
///
/// # Examples
///
/// ```
/// use molframe_core::contract::{Analysis, AnalysisPolicy, Coverage, Status};
///
/// let policy = AnalysisPolicy::default();
/// let radius = Analysis::complete(14.2_f64, Coverage::complete(1_960), &policy);
///
/// assert_eq!(*radius, 14.2);
/// assert_eq!(radius.status, Status::Complete);
/// assert_eq!(radius.coverage.fraction(), 1.0);
/// ```
#[derive(Clone, Debug)]
pub struct Analysis<T> {
    /// The answer.
    pub value: T,
    /// How far the analysis got.
    pub status: Status,
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
    /// A complete result.
    #[must_use]
    pub fn complete(value: T, coverage: Coverage, policy: &AnalysisPolicy) -> Self {
        Self {
            value,
            status: Status::Complete,
            coverage,
            warnings: Vec::new(),
            assumptions: Vec::new(),
            provenance: Provenance::new(policy),
        }
    }

    /// A result computed over less than it intended.
    #[must_use]
    pub fn partial(value: T, coverage: Coverage, policy: &AnalysisPolicy) -> Self {
        Self {
            status: Status::Partial,
            ..Self::complete(value, coverage, policy)
        }
    }

    /// A refusal to answer.
    ///
    /// The value is whatever stands in for "nothing"; the status is the result.
    #[must_use]
    pub fn indeterminate(value: T, coverage: Coverage, policy: &AnalysisPolicy) -> Self {
        Self {
            status: Status::Indeterminate,
            ..Self::complete(value, coverage, policy)
        }
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

    /// The decisions the caller never made that could have changed the answer.
    pub fn silent_assumptions(&self) -> impl Iterator<Item = &Assumption> {
        self.assumptions
            .iter()
            .filter(|assumption| assumption.is_silent_and_material())
    }

    /// Applies a function to the value, keeping everything else.
    pub fn map<U>(self, transform: impl FnOnce(T) -> U) -> Analysis<U> {
        Analysis {
            value: transform(self.value),
            status: self.status,
            coverage: self.coverage,
            warnings: self.warnings,
            assumptions: self.assumptions,
            provenance: self.provenance,
        }
    }
}

impl<T> Deref for Analysis<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.value
    }
}

#[cfg(test)]
#[path = "analysis_tests.rs"]
mod tests;
