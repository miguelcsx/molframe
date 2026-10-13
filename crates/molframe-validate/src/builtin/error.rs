//! Failure modes of the built-in validation.

use molframe_core::parallel::WorkerPanicked;
use molframe_spatial::SpatialError;

use molframe_chem::TableError;

/// Why a built-in validation could not complete.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum BuiltinError {
    /// A reference table is missing, unreadable or malformed.
    #[error(transparent)]
    Table(#[from] TableError),
    /// A worker thread panicked.
    #[error(transparent)]
    Worker(#[from] WorkerPanicked),
    /// The neighbour search behind the clash score failed.
    #[error("clash search failed: {0}")]
    Spatial(SpatialErrorKind),
}

/// A copyable summary of a [`SpatialError`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpatialErrorKind;

impl std::fmt::Display for SpatialErrorKind {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("spatial search error")
    }
}

impl From<SpatialError> for BuiltinError {
    fn from(_: SpatialError) -> Self {
        Self::Spatial(SpatialErrorKind)
    }
}
