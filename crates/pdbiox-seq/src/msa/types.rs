//! Public controls and failures for multiple-sequence alignment.

use crate::{AlignError, Scoring};

/// Explicit controls for progressive alignment and refinement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MsaOptions {
    /// Residue substitution and affine-gap scores.
    pub scoring: Scoring,
    /// Number of deterministic leave-one-out refinement passes.
    pub refinement_passes: usize,
    /// Hard ceiling for guide, dynamic-programming, and alignment storage.
    pub memory_limit_bytes: usize,
}

impl MsaOptions {
    /// Default operation-owned workspace ceiling.
    pub const DEFAULT_MEMORY_LIMIT_BYTES: usize = 100_000_000;

    /// Creates a single-linkage progressive alignment without refinement.
    #[must_use]
    pub const fn progressive(scoring: Scoring) -> Self {
        Self {
            scoring,
            refinement_passes: 0,
            memory_limit_bytes: Self::DEFAULT_MEMORY_LIMIT_BYTES,
        }
    }

    /// Enables a fixed number of leave-one-out refinement passes.
    #[must_use]
    pub const fn with_refinement_passes(mut self, passes: usize) -> Self {
        self.refinement_passes = passes;
        self
    }

    /// Sets the hard workspace and output-storage ceiling.
    #[must_use]
    pub const fn with_memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit_bytes = bytes;
        self
    }
}

/// Invalid input or scoring policy for multiple-sequence alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MsaError {
    /// Gap scores must be non-positive so opening or extending cannot reward a gap.
    InvalidGapScore,
    /// The requested ceiling is zero or exceeds the hard 500 MB boundary.
    InvalidMemoryLimit {
        /// Caller-selected ceiling.
        requested: usize,
    },
    /// The number of rows or columns exceeded representable indexing.
    DimensionOverflow,
    /// An exact pairwise or profile score exceeded its numeric domain.
    NumericOverflow,
    /// The requested guide or alignment cannot fit the explicit memory ceiling.
    MemoryLimit {
        /// Conservative bytes required at the rejected allocation boundary.
        required: usize,
        /// Caller-selected ceiling.
        limit: usize,
    },
}

impl std::fmt::Display for MsaError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidGapScore => {
                formatter.write_str("MSA affine-gap scores must be non-positive")
            }
            Self::InvalidMemoryLimit { requested } => write!(
                formatter,
                "MSA memory limit {requested} must be at least one byte"
            ),
            Self::DimensionOverflow => {
                formatter.write_str("MSA dimensions exceed the supported range")
            }
            Self::NumericOverflow => {
                formatter.write_str("MSA score exceeds the supported exact numeric range")
            }
            Self::MemoryLimit { required, limit } => {
                write!(
                    formatter,
                    "MSA requires {required} bytes, over the {limit} byte memory limit"
                )
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
