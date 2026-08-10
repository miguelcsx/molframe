//! Errors that preserve setup, frame and kernel failures.

use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_traj::TrajectoryError;
use std::fmt;

/// Failure while adapting a structure kernel to governed frame execution.
#[derive(Debug)]
#[non_exhaustive]
pub enum GovernedAnalysisError<E> {
    /// The trajectory executor rejected the requested execution.
    Trajectory(TrajectoryError),
    /// A frame could not form a valid immutable structure snapshot.
    InvalidStructure(Vec<Diagnostic>),
    /// The kernel failed and its typed error is preserved.
    Kernel(E),
    /// The policy requires failure when intended data are missing.
    MissingData {
        /// Number of intended inputs that were absent.
        missing: u32,
        /// Number of intended inputs that were ambiguous.
        ambiguous: u32,
    },
    /// Per-frame coverage cannot be represented by the public coverage counters.
    CoverageOverflow,
    /// A single-frame entry point received a policy selecting multiple models.
    MultipleModelsRequested,
    /// Deterministic execution did not produce the one requested frame value.
    MissingFrameOutput,
    /// This crate version does not understand a newer non-exhaustive policy value.
    UnsupportedPolicyValue(&'static str),
    /// A kernel reported mutually inconsistent coverage counters.
    InvalidCoverage {
        /// Inputs the kernel intended to use.
        intended: u32,
        /// Inputs the kernel reported as used.
        used: u32,
        /// Inputs the kernel reported as absent.
        missing: u32,
        /// Inputs the kernel reported as ambiguous.
        ambiguous: u32,
    },
}

impl<E> From<TrajectoryError> for GovernedAnalysisError<E> {
    fn from(error: TrajectoryError) -> Self {
        Self::Trajectory(error)
    }
}

impl<E: fmt::Display> fmt::Display for GovernedAnalysisError<E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Trajectory(error) => error.fmt(formatter),
            Self::InvalidStructure(findings) => {
                write!(formatter, "frame produced {} structure diagnostics", findings.len())
            }
            Self::Kernel(error) => write!(formatter, "analysis kernel failed: {error}"),
            Self::MissingData { missing, ambiguous } => write!(
                formatter,
                "analysis policy rejects {missing} missing and {ambiguous} ambiguous inputs"
            ),
            Self::CoverageOverflow => formatter.write_str("analysis coverage exceeds u32 limits"),
            Self::MultipleModelsRequested => formatter.write_str(
                "single-structure analysis requires model First or Index; use trajectory analysis for All or Ensemble",
            ),
            Self::MissingFrameOutput => {
                formatter.write_str("single-frame analysis produced no frame value")
            }
            Self::UnsupportedPolicyValue(field) => {
                write!(formatter, "unsupported analysis policy value for {field}")
            }
            Self::InvalidCoverage {
                intended,
                used,
                missing,
                ambiguous,
            } => write!(
                formatter,
                "invalid coverage: intended {intended}, used {used}, missing {missing}, ambiguous {ambiguous}"
            ),
        }
    }
}

impl<E> std::error::Error for GovernedAnalysisError<E> where E: std::error::Error + 'static {}
