//! Nearest neighbours of a separate query cloud among a target cloud.

use super::knn_within;
use crate::{KnnOptions, KnnResult, PeriodicBox, SpatialBackend, SpatialError};
use molframe_core::{AtomSelection, ExecutionContext};

/// Finds the `k` nearest `targets` of every point of `queries`.
///
/// Rows follow `queries`; neighbour indices refer to `targets`. Semantics
/// otherwise match [`knn_within`].
///
/// # Errors
///
/// Returns an error when the combined point count exceeds a 32-bit index, or
/// for any reason [`knn_within`] does.
#[allow(clippy::too_many_arguments)]
pub fn knn_between(
    targets: &[[f32; 3]],
    queries: &[[f32; 3]],
    k: usize,
    options: KnnOptions,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<KnnResult, SpatialError> {
    let total = targets
        .len()
        .checked_add(queries.len())
        .and_then(|total| u32::try_from(total).ok())
        .ok_or(SpatialError::AtomOutOfBounds(u32::MAX))?;
    let split =
        u32::try_from(targets.len()).map_err(|_| SpatialError::AtomOutOfBounds(u32::MAX))?;
    let mut joined = Vec::with_capacity(targets.len() + queries.len());
    joined.extend_from_slice(targets);
    joined.extend_from_slice(queries);
    knn_within(
        &joined,
        &AtomSelection::range(split..total),
        &AtomSelection::range(0..split),
        k,
        options,
        backend,
        periodic,
        context,
    )
}
