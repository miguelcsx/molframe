//! The words a spatial backend is spelled with.

use super::SpatialBackend;
use molframe_core::contract::{PolicyParseError, canonical_spelling};
use std::fmt;
use std::str::FromStr;

impl SpatialBackend {
    /// The words a backend can be spelled with, in declaration order.
    pub const NAMES: &'static [&'static str] = &[
        "brute-force",
        "cell-list",
        "kd-tree",
        "neighbor-list",
        "auto",
    ];

    /// The canonical spelling of this backend.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::BruteForce => "brute-force",
            Self::CellList => "cell-list",
            Self::KdTree => "kd-tree",
            Self::NeighborList => "neighbor-list",
            Self::Auto => "auto",
        }
    }
}

impl fmt::Display for SpatialBackend {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl FromStr for SpatialBackend {
    type Err = PolicyParseError;

    /// `cell` is accepted for `cell-list`, and an underscore for a hyphen.
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match canonical_spelling(value).as_str() {
            "brute-force" => Ok(Self::BruteForce),
            "cell-list" | "cell" => Ok(Self::CellList),
            "kd-tree" => Ok(Self::KdTree),
            "neighbor-list" => Ok(Self::NeighborList),
            "auto" => Ok(Self::Auto),
            _ => Err(PolicyParseError::new(
                "backend",
                value,
                &Self::NAMES.join(", "),
            )),
        }
    }
}
