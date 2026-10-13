//! Named input sizes, so a report column means the same thing everywhere.

use std::fmt;

/// A synthetic input size class, by atom count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Tier {
    /// One thousand atoms.
    Tiny,
    /// Ten thousand atoms: a large single chain.
    Small,
    /// One hundred thousand atoms: a large complex.
    Medium,
    /// One million atoms: a viral capsid or a solvated box.
    Large,
    /// Ten million atoms.
    Huge,
}

/// An unrecognised tier name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TierParseError(pub String);

impl fmt::Display for TierParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "unknown tier `{}`", self.0)
    }
}

impl std::error::Error for TierParseError {}

impl Tier {
    /// Every tier, smallest first.
    pub const ALL: [Tier; 5] = [
        Tier::Tiny,
        Tier::Small,
        Tier::Medium,
        Tier::Large,
        Tier::Huge,
    ];

    /// The number of atoms this tier stands for.
    #[must_use]
    pub const fn atoms(self) -> usize {
        match self {
            Tier::Tiny => 1_000,
            Tier::Small => 10_000,
            Tier::Medium => 100_000,
            Tier::Large => 1_000_000,
            Tier::Huge => 10_000_000,
        }
    }

    /// The stable lower-case name used on the command line and in reports.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Tier::Tiny => "tiny",
            Tier::Small => "small",
            Tier::Medium => "medium",
            Tier::Large => "large",
            Tier::Huge => "huge",
        }
    }

    /// Parses a label produced by [`Tier::label`].
    ///
    /// # Errors
    ///
    /// Returns [`TierParseError`] for any other text.
    pub fn parse(text: &str) -> Result<Self, TierParseError> {
        Self::ALL
            .into_iter()
            .find(|tier| tier.label() == text)
            .ok_or_else(|| TierParseError(text.to_owned()))
    }
}
