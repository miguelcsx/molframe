//! Public spatial query values.

use pdbiox_core::{Code, Diagnostic};
use std::fmt;

/// A deterministic spatial implementation choice.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[non_exhaustive]
pub enum SpatialBackend {
    /// Direct comparison, for small one-shot workloads.
    BruteForce,
    /// Contiguous fixed-radius grid.
    CellList,
    /// Balanced static k-dimensional tree.
    KdTree,
    /// Verlet candidate list with a displacement skin.
    NeighborList,
    /// Select from workload shape.
    Auto,
}

/// One unique unordered atom pair within a cutoff.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct NeighborPair {
    /// Lower atom index.
    pub first: u32,
    /// Higher atom index.
    pub second: u32,
    /// Squared separation in ångström squared.
    pub distance_squared: f32,
}

impl NeighborPair {
    pub(crate) fn new(left: u32, right: u32, distance_squared: f32) -> Self {
        Self {
            first: left.min(right),
            second: left.max(right),
            distance_squared,
        }
    }
}

/// Why a spatial query could not be evaluated.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum SpatialError {
    /// The cutoff was negative or non-finite.
    InvalidCutoff,
    /// A selection named an atom outside the coordinate array.
    AtomOutOfBounds(u32),
    /// Cell parameters do not define an invertible three-dimensional box.
    InvalidCell,
    /// Positions moved beyond the displacement skin used to build a list.
    StaleNeighborList,
}

impl SpatialError {
    /// Converts the backend error into pdbiox's registered diagnostic model.
    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        match self {
            Self::InvalidCutoff => Diagnostic::new(Code::E4002),
            Self::AtomOutOfBounds(atom) => {
                Diagnostic::new(Code::E6009).with_context("atom", atom.to_string())
            }
            Self::InvalidCell => Diagnostic::new(Code::E5004),
            Self::StaleNeighborList => Diagnostic::new(Code::E9001),
        }
    }
}

impl fmt::Display for SpatialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCutoff => f.write_str("cutoff must be finite and non-negative"),
            Self::AtomOutOfBounds(atom) => write!(f, "atom {atom} is outside the coordinate array"),
            Self::InvalidCell => f.write_str("unit cell is degenerate or non-finite"),
            Self::StaleNeighborList => {
                f.write_str("positions moved beyond the neighbour-list displacement skin")
            }
        }
    }
}

impl std::error::Error for SpatialError {}
