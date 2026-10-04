//! Why a policy space could not be planned.

use molframe_core::contract::PolicyField;
use std::fmt;

/// Why a policy space could not be planned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum PlanError {
    /// A named dimension has no alternatives.
    EmptyDimension(PolicyField),
    /// A field was declared more than once.
    DuplicateDimension(PolicyField),
    /// The Cartesian product overflowed the platform's index size.
    CostOverflow,
    /// A run of the space combines decisions that contradict each other.
    ///
    /// Restrict the dimensions so every combination is a system that can run;
    /// the policy at that index says which pair. `AnalysisPolicy::check_consistency`
    /// names it.
    Conflict {
        /// Zero-based run whose policy is contradictory.
        run: usize,
    },
    /// A varied field is one the analysis never applied, so varying it measures nothing.
    NotRead(PolicyField),
    /// Every combination was removed as contradictory or forbidden.
    NoUniverse,
    /// The requested space exceeds the caller's bound.
    LimitExceeded {
        /// Exact requested number of runs.
        cost: usize,
        /// Configured upper bound.
        limit: usize,
    },
}

impl fmt::Display for PlanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDimension(field) => {
                write!(formatter, "{} has no alternatives", field.name())
            }
            Self::DuplicateDimension(field) => {
                write!(formatter, "{} is varied twice", field.name())
            }
            Self::CostOverflow => formatter.write_str("policy-space cost overflowed usize"),
            Self::NoUniverse => {
                formatter.write_str("every combination is contradictory or forbidden")
            }
            Self::Conflict { run } => write!(
                formatter,
                "run {run} combines decisions that contradict each other"
            ),
            Self::NotRead(field) => write!(
                formatter,
                "the analysis never reads {}, so varying it would report stability it has not earned",
                field.name()
            ),
            Self::LimitExceeded { cost, limit } => {
                write!(
                    formatter,
                    "policy space needs {cost} runs, exceeding limit {limit}"
                )
            }
        }
    }
}

impl std::error::Error for PlanError {}
