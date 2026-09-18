//! Stable diagnostic codes and the two axes that govern them.
//!
//! Codes are stable forever. Users grep for them and scripts branch on them, so
//! a withdrawn code is retired rather than reused. The number's leading digit
//! names the class of problem, which lets a caller filter for, say, every
//! conversion loss without enumerating individual codes.
//!
//! Severity and strictness are separate axes on purpose. Severity says how bad a
//! finding is; strictness says how the caller wants to react to it. Keeping them
//! apart is what lets one implementation serve a curation pipeline that must
//! reject anything irregular and an analysis script that has to cope with what
//! the archive actually contains.

use std::fmt;

/// Whether a code denotes a failure or an observation.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Kind {
    /// An error: rendered `MOLFRAME-Ennnn`.
    Error,
    /// A warning: rendered `MOLFRAME-Wnnnn`.
    Warning,
}

impl Kind {
    const fn letter(self) -> char {
        match self {
            Self::Error => 'E',
            Self::Warning => 'W',
        }
    }
}

/// How bad a finding is, independent of how the caller wants to react.
///
/// Ordered so that greater is worse, which makes "at least this bad" a
/// comparison and makes sorting diagnostics worst-first a plain sort.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Severity {
    /// Worth recording, and nothing more.
    Info,
    /// Unusual, but probably what the depositor intended.
    Loose,
    /// Violates the specification; usually still interpretable.
    Strict,
    /// The result computed from this would be wrong.
    Invalidating,
    /// The data cannot be interpreted at all.
    Breaking,
}

/// How much irregularity the caller is willing to accept.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum Strictness {
    /// Anything violating the specification is an error.
    Strict,
    /// Anything that would make the result wrong is an error.
    #[default]
    Medium,
    /// Only uninterpretable data is an error.
    Loose,
}

impl Severity {
    /// Returns true when a finding of this severity is an error under
    /// `strictness`.
    ///
    /// # Examples
    ///
    /// ```
    /// use molframe_core::{Severity, Strictness};
    ///
    /// assert!(Severity::Strict.is_error(Strictness::Strict));
    /// assert!(!Severity::Strict.is_error(Strictness::Medium));
    /// assert!(Severity::Breaking.is_error(Strictness::Loose));
    /// ```
    #[must_use]
    pub const fn is_error(self, strictness: Strictness) -> bool {
        let floor = match strictness {
            Strictness::Strict => Self::Strict,
            Strictness::Medium => Self::Invalidating,
            Strictness::Loose => Self::Breaking,
        };
        (self as u8) >= (floor as u8)
    }

    /// The word used when rendering a finding of this severity.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Info => "note",
            Self::Loose | Self::Strict => "warning",
            Self::Invalidating | Self::Breaking => "error",
        }
    }
}

/// A stable diagnostic code, rendered as `MOLFRAME-E1103` or `MOLFRAME-W3011`.
///
/// # Examples
///
/// ```
/// use molframe_core::Code;
///
/// assert_eq!(Code::E1103.to_string(), "MOLFRAME-E1103");
/// assert_eq!(Code::E1103.cause(), "loop_ row count is not a multiple of the column count");
/// assert!(!Code::E1103.remedy().is_empty());
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Code {
    kind: Kind,
    number: u16,
}

impl Code {
    /// Creates a code from its kind and number.
    ///
    /// Prefer the associated constants; this exists so the registry can be
    /// declared as data.
    #[must_use]
    pub const fn new(kind: Kind, number: u16) -> Self {
        Self { kind, number }
    }

    /// Whether this code is an error or a warning.
    #[must_use]
    pub const fn kind(self) -> Kind {
        self.kind
    }

    /// The numeric part of the code.
    #[must_use]
    pub const fn number(self) -> u16 {
        self.number
    }

    /// The class this code belongs to, taken from its leading digit.
    #[must_use]
    pub const fn class(self) -> Class {
        Class::of(self.number)
    }

    fn entry(self) -> Option<&'static Entry> {
        let registry = super::registry::ENTRIES;
        let found = registry
            .binary_search_by(|entry| entry.code.cmp(&self))
            .ok()?;
        registry.get(found)
    }

    /// The registered severity, or [`Severity::Breaking`] for a code that is not
    /// in the registry.
    ///
    /// An unregistered code can only arise from a programming mistake, so the
    /// conservative reading is the right one.
    #[must_use]
    pub fn severity(self) -> Severity {
        self.entry()
            .map_or(Severity::Breaking, |entry| entry.severity)
    }

    /// What went wrong.
    #[must_use]
    pub fn cause(self) -> &'static str {
        self.entry()
            .map_or("unregistered diagnostic code", |entry| entry.cause)
    }

    /// What to do about it.
    ///
    /// Every registered code has one. A diagnostic a user cannot act on is an
    /// incomplete diagnostic, so the remedy is part of the registry rather than
    /// something each call site is trusted to remember.
    #[must_use]
    pub fn remedy(self) -> &'static str {
        self.entry()
            .map_or("report this as a molframe bug", |entry| entry.remedy)
    }

    /// Returns true when this code is present in the registry.
    #[must_use]
    pub fn is_registered(self) -> bool {
        self.entry().is_some()
    }

    /// Every registered code, ordered by kind and then number.
    pub fn registered() -> impl Iterator<Item = Self> {
        super::registry::ENTRIES.iter().map(|entry| entry.code)
    }
}

impl fmt::Display for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MOLFRAME-{}{:04}", self.kind.letter(), self.number)
    }
}

impl fmt::Debug for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self}")
    }
}

/// The class of problem a code's leading digit names.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[non_exhaustive]
pub enum Class {
    /// Syntax and lexical structure.
    Syntax,
    /// Schema and dictionary conformance.
    Schema,
    /// Structural consistency of the model.
    Consistency,
    /// Conversion and data loss.
    Conversion,
    /// Geometry and numerics.
    Geometry,
    /// Policy and analysis contracts.
    Policy,
    /// Input, output and resources.
    Resource,
    /// An internal invariant was violated. These are bugs in molframe.
    Internal,
}

impl Class {
    const fn of(number: u16) -> Self {
        match number / 1000 {
            1 => Self::Syntax,
            2 => Self::Schema,
            3 => Self::Consistency,
            4 => Self::Conversion,
            5 => Self::Geometry,
            6 => Self::Policy,
            7 => Self::Resource,
            _ => Self::Internal,
        }
    }
}

/// One row of the code registry.
pub(super) struct Entry {
    pub(super) code: Code,
    pub(super) severity: Severity,
    pub(super) cause: &'static str,
    pub(super) remedy: &'static str,
}

#[cfg(test)]
#[path = "code_tests.rs"]
mod tests;
