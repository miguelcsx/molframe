//! Public spatial query values.

use molframe_core::{Code, Diagnostic};
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
    /// Shared execution memory cannot admit the complete spatial workspace.
    Memory(molframe_core::MemoryBudgetError),
    /// Cooperative cancellation stopped the spatial operation.
    Cancelled,
    /// The cutoff was negative or non-finite.
    InvalidCutoff,
    /// A selection named an atom outside the coordinate array.
    AtomOutOfBounds(u32),
    /// Cell parameters do not define an invertible three-dimensional box.
    InvalidCell,
    /// Positions moved beyond the displacement skin used to build a list.
    StaleNeighborList,
    /// A public planning option is invalid.
    InvalidOption(SpatialOption),
    /// A periodic k-d query exceeds its explicit image budget.
    PeriodicImageLimitExceeded {
        /// Images required for an exact query.
        required: usize,
        /// Maximum images permitted by the profile.
        maximum: usize,
    },
    /// A derived coordinate or index cannot be represented by the public type.
    NumericRangeExceeded,
    /// A scoped worker thread panicked during a parallel search.
    WorkerPanicked,
}

/// Identifies an invalid field in a spatial planning profile.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[non_exhaustive]
pub enum SpatialOption {
    /// Direct-comparison workload threshold.
    BruteForcePairLimit,
    /// Minimum target count for a k-d tree.
    KdTargetMinimum,
    /// Target-to-query ratio for a k-d tree.
    KdQueryRatio,
    /// Concrete periodic backend used by automatic planning.
    PeriodicBackend,
    /// Cutoff fraction used to derive a neighbour skin.
    NeighborSkinRatio,
    /// Minimum or absolute neighbour skin.
    NeighborSkinMinimum,
    /// Maximum allocated cell count.
    MaximumCellCount,
    /// Cell-edge coarsening multiplier.
    CellGrowthFactor,
    /// Maximum periodic k-d query image count.
    KdPeriodicImageLimit,
}

impl SpatialError {
    /// Converts the backend error into molframe's registered diagnostic model.
    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        match self {
            Self::Memory(error) => {
                Diagnostic::new(Code::E1901).with_context("execution_memory", error.to_string())
            }
            Self::Cancelled => {
                Diagnostic::new(Code::E1904).with_context("spatial_operation", "cancelled")
            }
            Self::InvalidCutoff => Diagnostic::new(Code::E4002),
            Self::AtomOutOfBounds(atom) => {
                Diagnostic::new(Code::E6009).with_context("atom", atom.to_string())
            }
            Self::InvalidCell => Diagnostic::new(Code::E5004),
            Self::StaleNeighborList => Diagnostic::new(Code::E9001),
            Self::InvalidOption(option) => {
                Diagnostic::new(Code::E4002).with_context("spatial_option", option.to_string())
            }
            Self::PeriodicImageLimitExceeded { required, maximum } => Diagnostic::new(Code::E4002)
                .with_context("required_periodic_images", required.to_string())
                .with_context("maximum_periodic_images", maximum.to_string()),
            Self::NumericRangeExceeded => {
                Diagnostic::new(Code::E4002).with_context("spatial_numeric_range", "exceeded")
            }
            Self::WorkerPanicked => {
                Diagnostic::new(Code::E4002).with_context("spatial_worker", "panicked")
            }
        }
    }
}

impl fmt::Display for SpatialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Memory(error) => error.fmt(f),
            Self::Cancelled => f.write_str("spatial operation cancelled"),
            Self::InvalidCutoff => f.write_str("cutoff must be finite and non-negative"),
            Self::AtomOutOfBounds(atom) => write!(f, "atom {atom} is outside the coordinate array"),
            Self::InvalidCell => f.write_str("unit cell is degenerate or non-finite"),
            Self::StaleNeighborList => {
                f.write_str("positions moved beyond the neighbour-list displacement skin")
            }
            Self::InvalidOption(option) => write!(f, "invalid spatial option: {option}"),
            Self::PeriodicImageLimitExceeded { required, maximum } => write!(
                f,
                "periodic k-d query requires {required} images but the configured limit is {maximum}"
            ),
            Self::NumericRangeExceeded => {
                f.write_str("a derived spatial value exceeds its representable numeric range")
            }
            Self::WorkerPanicked => f.write_str("a spatial worker thread panicked"),
        }
    }
}

impl std::error::Error for SpatialError {}

impl fmt::Display for SpatialOption {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::BruteForcePairLimit => "brute_force_pair_limit",
            Self::KdTargetMinimum => "kd_target_minimum",
            Self::KdQueryRatio => "kd_query_ratio",
            Self::PeriodicBackend => "periodic_backend",
            Self::NeighborSkinRatio => "neighbor_skin.cutoff_ratio",
            Self::NeighborSkinMinimum => "neighbor_skin.minimum",
            Self::MaximumCellCount => "cell_grid.maximum_cell_count",
            Self::CellGrowthFactor => "cell_grid.edge_growth_factor",
            Self::KdPeriodicImageLimit => "kd_periodic.maximum_image_count",
        })
    }
}

impl From<molframe_core::MemoryBudgetError> for SpatialError {
    fn from(error: molframe_core::MemoryBudgetError) -> Self {
        Self::Memory(error)
    }
}

#[cfg(test)]
#[path = "types_tests.rs"]
mod tests;
