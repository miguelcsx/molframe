//! Comparing two snapshots.

use molframe_core::structure::{
    DifferenceError, StructureDifference, StructureDifferenceOptions,
    structure_difference as engine_difference,
};

use super::handle::Structure;

/// Reports how two snapshots differ.
///
/// The comparison is over the topology tables and the first model's
/// coordinates; `coordinate_tolerance` is the absolute displacement below which
/// two positions count as the same.
///
/// # Errors
///
/// Returns [`DifferenceError::InvalidCoordinateTolerance`] when the tolerance is
/// negative or not finite.
pub fn structure_difference(
    left: &Structure,
    right: &Structure,
    options: StructureDifferenceOptions,
) -> Result<StructureDifference, DifferenceError> {
    engine_difference(left.engine(), right.engine(), options)
}
