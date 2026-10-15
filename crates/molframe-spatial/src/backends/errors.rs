//! Validation and error mapping shared by the index backends.

use crate::SpatialError;
use molframe_core::parallel::BlockExecutionError;

/// Validates a non-negative finite cutoff.
///
/// Runtime and auxiliary space are `O(1)`.
pub(super) fn validate_cutoff(cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

/// Validates atom indices against a coordinate-array length.
///
/// Runtime is `O(N)` and no allocation is performed.
pub(super) fn validate_indices(atoms: &[u32], position_count: usize) -> Result<(), SpatialError> {
    for &atom in atoms {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if index >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom));
        }
    }

    Ok(())
}

/// Maps a block-execution failure onto the spatial error vocabulary.
pub(super) fn block_error(error: &BlockExecutionError<SpatialError>) -> SpatialError {
    match error {
        BlockExecutionError::Memory(error) => SpatialError::Memory(*error),
        BlockExecutionError::Cancelled => SpatialError::Cancelled,
        BlockExecutionError::Worker(_) => SpatialError::WorkerPanicked,
        BlockExecutionError::Operation(error) => *error,
    }
}
