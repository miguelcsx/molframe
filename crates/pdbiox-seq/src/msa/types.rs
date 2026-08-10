//! Public controls and failures for multiple-sequence alignment.

use crate::{AlignError, Scoring};

/// Explicit controls for progressive alignment and refinement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MsaOptions {
    /// Residue substitution and affine-gap scores.
    pub scoring: Scoring,
    /// Number of deterministic leave-one-out refinement passes.
    pub refinement_passes: usize,
}

impl MsaOptions {
    /// Creates a progressive alignment without refinement.
    #[must_use]
    pub const fn progressive(scoring: Scoring) -> Self {
        Self {
            scoring,
            refinement_passes: 0,
        }
    }

    /// Enables a fixed number of leave-one-out refinement passes.
    #[must_use]
    pub const fn with_refinement_passes(mut self, passes: usize) -> Self {
        self.refinement_passes = passes;
        self
    }
}

/// Invalid input or scoring policy for multiple-sequence alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MsaError {
    /// Gap scores must be non-positive so opening or extending cannot reward a gap.
    InvalidGapScore,
    /// The number of rows or columns exceeded representable indexing.
    DimensionOverflow,
    /// An exact pairwise or profile score exceeded its numeric domain.
    NumericOverflow,
}

impl std::fmt::Display for MsaError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidGapScore => {
                formatter.write_str("MSA affine-gap scores must be non-positive")
            }
            Self::DimensionOverflow => {
                formatter.write_str("MSA dimensions exceed the supported range")
            }
            Self::NumericOverflow => {
                formatter.write_str("MSA score exceeds the supported exact numeric range")
            }
        }
    }
}

impl std::error::Error for MsaError {}

impl From<AlignError> for MsaError {
    fn from(_: AlignError) -> Self {
        Self::NumericOverflow
    }
}
