//! A finished set of findings, as the error an operation returns.

use super::code::Code;
use super::finding::Diagnostic;
use std::fmt;

/// The findings an operation produced, as a finished set.
///
/// [`Diagnostics`](super::Diagnostics) accumulates and
/// [`Diagnostics::finish`](super::Diagnostics::finish) imposes the order; this is
/// what comes out the other side, and it is an error because reaching it means the
/// operation did not produce what was asked for. The findings say why, in the
/// order two runs agree on.
///
/// It derefs to the findings, so a caller that wants the list takes it as a
/// slice; `for finding in &set` iterates without naming [`Deref`] at all.
///
/// # Examples
///
/// ```
/// use molframe_core::{Code, Diagnostic, Findings};
///
/// let findings = Findings::from(Diagnostic::new(Code::E1001));
/// assert_eq!(findings.len(), 1);
/// assert_eq!(findings.code(), Some(Code::E1001));
///
/// // The order `Diagnostics::finish` produces, and a report that reads it.
/// println!("{findings}");
/// ```
///
/// [`Deref`]: std::ops::Deref
#[derive(Clone, Debug)]
pub struct Findings(Vec<Diagnostic>);

impl Findings {
    /// The code of the finding that stopped the operation, if any.
    ///
    /// A stable identity to branch on, for a caller that would rather not match
    /// on message text. `None` when the set is empty — a set that stopped
    /// nothing rather than a set that failed to say why.
    #[must_use]
    pub fn code(&self) -> Option<Code> {
        self.0.first().map(Diagnostic::code)
    }

    /// The findings.
    ///
    /// Declared rather than inherited from [`Deref`], whose `[Diagnostic]` target
    /// reaches `[T]::as_slice` — an unstable inherent method that shadows a stable
    /// reading of the same call. [`Diagnostics`](super::Diagnostics) names it the
    /// same way.
    ///
    /// [`Deref`]: std::ops::Deref
    #[must_use]
    pub fn as_slice(&self) -> &[Diagnostic] {
        &self.0
    }

    /// Consumes the set, returning the findings.
    #[must_use]
    pub fn into_vec(self) -> Vec<Diagnostic> {
        self.0
    }
}

impl std::ops::Deref for Findings {
    type Target = [Diagnostic];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> IntoIterator for &'a Findings {
    type Item = &'a Diagnostic;
    type IntoIter = std::slice::Iter<'a, Diagnostic>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

/// Prints the first finding and how many follow it.
///
/// A set can hold thousands, and an error message that renders all of them is
/// one nobody reads. The first is the worst by
/// [`Diagnostics::finish`](super::Diagnostics::finish)'s order, which is the one
/// worth naming.
impl fmt::Display for Findings {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0.as_slice() {
            [] => f.write_str("no findings"),
            [only] => write!(f, "{only}"),
            [first, rest @ ..] => write!(f, "{first} (and {} more)", rest.len()),
        }
    }
}

impl std::error::Error for Findings {}

impl From<Diagnostic> for Findings {
    fn from(finding: Diagnostic) -> Self {
        Self(vec![finding])
    }
}

impl From<Vec<Diagnostic>> for Findings {
    fn from(findings: Vec<Diagnostic>) -> Self {
        Self(findings)
    }
}

impl From<Findings> for Vec<Diagnostic> {
    fn from(findings: Findings) -> Self {
        findings.0
    }
}

#[cfg(test)]
#[path = "findings_tests.rs"]
mod tests;
