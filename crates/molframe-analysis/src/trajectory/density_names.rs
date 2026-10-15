//! The words a Cartesian axis is spelled with.

use crate::CartesianAxis;
use molframe_core::contract::{PolicyParseError, canonical_spelling};
use std::fmt;
use std::str::FromStr;

impl CartesianAxis {
    /// The words an axis can be spelled with, in declaration order.
    pub const NAMES: &'static [&'static str] = &["x", "y", "z"];

    /// The canonical spelling of this axis.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
        }
    }
}

impl fmt::Display for CartesianAxis {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl FromStr for CartesianAxis {
    type Err = PolicyParseError;

    /// Case-insensitive: `Z` is `z`.
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match canonical_spelling(value).to_ascii_lowercase().as_str() {
            "x" => Ok(Self::X),
            "y" => Ok(Self::Y),
            "z" => Ok(Self::Z),
            _ => Err(PolicyParseError::new(
                "axis",
                value,
                &Self::NAMES.join(", "),
            )),
        }
    }
}
