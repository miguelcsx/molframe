//! Strict TNG read and write failures.

/// Malformed, unsupported or unrepresentable TNG data.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum TngError {
    /// The native TNG implementation rejected the container or operation.
    #[error("TNG container error: {0}")]
    Container(String),
    /// A TNG invariant failed inside the native implementation.
    #[error("TNG container violated an internal invariant")]
    InternalInvariant,
    /// No position frames, no particles or contradictory shapes were found.
    #[error("invalid TNG trajectory shape")]
    InvalidShape,
    /// A coordinate, time, cell or unit conversion is not finite or representable.
    #[error("invalid TNG numeric value")]
    InvalidValue,
    /// Simulation steps must be non-negative and strictly increasing.
    #[error("TNG simulation steps are not uniformly increasing")]
    InvalidSteps,
    /// Frame times do not define one positive time interval per simulation step.
    #[error("TNG frame times are not uniformly increasing with simulation steps")]
    InvalidTimeAxis,
    /// An auxiliary array disagrees with the position atom count.
    #[error("TNG auxiliary array has {found} atoms; expected {expected}")]
    AtomCountMismatch {
        /// Position atom count.
        expected: usize,
        /// Auxiliary atom count.
        found: usize,
    },
}

impl From<tng_rs::TngError> for TngError {
    fn from(error: tng_rs::TngError) -> Self {
        Self::Container(error.to_string())
    }
}
