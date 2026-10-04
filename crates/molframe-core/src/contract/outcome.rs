//! The answer, or the reason there is none.
//!
//! A refusal to answer used to be a [`Status`](super::Status) beside a value
//! that "is not meaningful". A value that must not be read but can be is a
//! documented hazard, not a guarantee. Here the two states are different
//! shapes: [`Outcome::Determinate`] holds the value and [`Outcome::Indeterminate`]
//! holds why there is none, so there is nothing to misread.

use super::analysis::{Coverage, Status};
use super::policy::MissingPolicy;
use std::fmt;

/// Why no defensible answer exists under a policy.
#[derive(Clone, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Indeterminacy {
    /// The policy refuses to compute over atoms that are missing or ambiguous.
    MissingInputs {
        /// Intended atoms that were absent.
        missing: u32,
        /// Intended atoms that were present but not uniquely resolvable.
        ambiguous: u32,
    },
    /// The alternate-conformation rule could not select a consistent set of atoms.
    UnresolvedConformations,
    /// One frame of a series had no defensible answer, so the series has none.
    Frame {
        /// Zero-based frame index.
        frame: usize,
        /// Why that frame has none.
        reason: Box<Indeterminacy>,
    },
    /// A reason supplied by the analysis itself.
    Other(Box<str>),
}

impl fmt::Display for Indeterminacy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingInputs { missing, ambiguous } => write!(
                formatter,
                "the policy rejects {missing} missing and {ambiguous} ambiguous inputs"
            ),
            Self::UnresolvedConformations => {
                formatter.write_str("the alternate-conformation rule selected no consistent atoms")
            }
            Self::Frame { frame, reason } => write!(formatter, "frame {frame}: {reason}"),
            Self::Other(reason) => formatter.write_str(reason),
        }
    }
}

impl std::error::Error for Indeterminacy {}

/// A value, or the reason there is none.
///
/// There is deliberately no `Deref` and no accessor that returns a value for the
/// indeterminate case: the way to a number goes through a match, [`Outcome::value`]
/// or [`Outcome::into_result`], each of which says what happens otherwise.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Outcome<T> {
    /// A defensible answer.
    Determinate(T),
    /// No defensible answer exists.
    Indeterminate(Indeterminacy),
}

impl<T> Outcome<T> {
    /// Returns true when there is an answer.
    #[must_use]
    pub const fn is_determinate(&self) -> bool {
        matches!(self, Self::Determinate(_))
    }

    /// The answer, if there is one.
    #[must_use]
    pub const fn value(&self) -> Option<&T> {
        match self {
            Self::Determinate(value) => Some(value),
            Self::Indeterminate(_) => None,
        }
    }

    /// The reason there is no answer, if there is none.
    #[must_use]
    pub const fn indeterminacy(&self) -> Option<&Indeterminacy> {
        match self {
            Self::Determinate(_) => None,
            Self::Indeterminate(reason) => Some(reason),
        }
    }

    /// The answer, or the reason it is missing.
    ///
    /// # Errors
    ///
    /// Returns the [`Indeterminacy`] when no defensible answer exists.
    pub fn into_result(self) -> Result<T, Indeterminacy> {
        match self {
            Self::Determinate(value) => Ok(value),
            Self::Indeterminate(reason) => Err(reason),
        }
    }

    /// Borrows the answer, or the reason it is missing.
    ///
    /// # Errors
    ///
    /// Returns the [`Indeterminacy`] when no defensible answer exists.
    pub const fn as_result(&self) -> Result<&T, &Indeterminacy> {
        match self {
            Self::Determinate(value) => Ok(value),
            Self::Indeterminate(reason) => Err(reason),
        }
    }

    /// Transforms the answer and leaves a refusal as it is.
    pub fn map<U>(self, transform: impl FnOnce(T) -> U) -> Outcome<U> {
        match self {
            Self::Determinate(value) => Outcome::Determinate(transform(value)),
            Self::Indeterminate(reason) => Outcome::Indeterminate(reason),
        }
    }
}

/// How complete a determinate answer is.
///
/// The qualities of an answer that exists. "No answer" is not one of them: it is
/// [`Outcome::Indeterminate`].
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum Quality {
    /// Every intended atom was available and unambiguous.
    #[default]
    Complete,
    /// Computed, but some intended atoms were not available.
    Partial,
    /// Several defensible answers exist; one was chosen deterministically.
    Ambiguous,
}

impl Quality {
    /// The weaker of two qualities, as a combined result is only as good as its parts.
    #[must_use]
    pub fn worst(self, other: Self) -> Self {
        // Ambiguity is the stronger caveat: a partial answer is still unique.
        match (self, other) {
            (Self::Ambiguous, _) | (_, Self::Ambiguous) => Self::Ambiguous,
            (Self::Partial, _) | (_, Self::Partial) => Self::Partial,
            (Self::Complete, Self::Complete) => Self::Complete,
        }
    }
}

impl From<Quality> for Status {
    fn from(quality: Quality) -> Self {
        match quality {
            Quality::Complete => Self::Complete,
            Quality::Partial => Self::Partial,
            Quality::Ambiguous => Self::Ambiguous,
        }
    }
}

/// Why a missing-input policy could not be applied.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum MissingPolicyError {
    /// The policy is to fail when intended inputs are missing or ambiguous.
    Fail {
        /// Intended inputs that were absent.
        missing: u32,
        /// Intended inputs that were present but not uniquely resolvable.
        ambiguous: u32,
    },
}

impl fmt::Display for MissingPolicyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fail { missing, ambiguous } => write!(
                formatter,
                "the policy fails on {missing} missing and {ambiguous} ambiguous inputs"
            ),
        }
    }
}

impl std::error::Error for MissingPolicyError {}

/// What the policy for missing inputs makes of a value computed over `coverage`.
///
/// Complete coverage keeps the value at the quality it came with. Otherwise
/// `Ignore` keeps it as it is, `Report` marks it partial, `Indeterminate` drops
/// the value for a refusal that says what was missing, and `Fail` is an error.
/// This is the one place the rule lives; every analysis that governs its
/// completeness goes through it.
///
/// # Errors
///
/// Returns [`MissingPolicyError::Fail`] when the policy is to fail and inputs are
/// missing or ambiguous.
pub fn resolve_missing<T>(
    value: T,
    quality: Quality,
    coverage: Coverage,
    policy: MissingPolicy,
) -> Result<(Outcome<T>, Quality), MissingPolicyError> {
    if coverage.missing == 0 && coverage.ambiguous == 0 {
        return Ok((Outcome::Determinate(value), quality));
    }
    match policy {
        MissingPolicy::Ignore => Ok((Outcome::Determinate(value), quality)),
        MissingPolicy::Report => Ok((Outcome::Determinate(value), quality.worst(Quality::Partial))),
        MissingPolicy::Indeterminate => Ok((
            Outcome::Indeterminate(Indeterminacy::MissingInputs {
                missing: coverage.missing,
                ambiguous: coverage.ambiguous,
            }),
            Quality::Complete,
        )),
        MissingPolicy::Fail => Err(MissingPolicyError::Fail {
            missing: coverage.missing,
            ambiguous: coverage.ambiguous,
        }),
    }
}
