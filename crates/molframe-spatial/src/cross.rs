//! Neighbours between coordinate sets with independent index spaces.

use molframe_core::{AtomSelection, ExecutionContext};

use crate::{SpatialBackend, SpatialError, pairs_within};

/// A neighbour relation with indices local to each input coordinate set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrossPair {
    /// Index in the first coordinate set.
    pub first: u32,
    /// Index in the second coordinate set.
    pub second: u32,
    /// Squared separation in ångström squared.
    pub distance_squared: f32,
}

/// Find cross-set neighbours without producing within-set pairs.
///
/// The coordinate workspace respects the shared execution memory budget. Output
/// order follows the first index and then the second index deterministically.
///
/// # Errors
/// Returns an error for invalid coordinates or cutoff, overflowing indices,
/// cancellation, or insufficient execution memory.
pub fn cross_pairs(
    first: &[[f32; 3]],
    second: &[[f32; 3]],
    cutoff: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<Vec<CrossPair>, SpatialError> {
    let left_count = u32::try_from(first.len()).map_err(|_| SpatialError::NumericRangeExceeded)?;
    let count = first
        .len()
        .checked_add(second.len())
        .ok_or(SpatialError::NumericRangeExceeded)?;
    let total = u32::try_from(count).map_err(|_| SpatialError::NumericRangeExceeded)?;
    let bytes = count
        .checked_mul(size_of::<[f32; 3]>())
        .ok_or(SpatialError::NumericRangeExceeded)?;
    let _reservation = context.try_reserve(bytes)?;
    let mut positions = Vec::with_capacity(count);
    positions.extend_from_slice(first);
    positions.extend_from_slice(second);
    let pairs = pairs_within(
        &positions,
        &AtomSelection::range(0..left_count),
        &AtomSelection::range(left_count..total),
        cutoff,
        backend,
        None,
        context,
    )?;
    Ok(pairs
        .into_iter()
        .map(|pair| CrossPair {
            first: pair.first,
            second: pair.second - left_count,
            distance_squared: pair.distance_squared,
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_indices_remain_local_even_when_the_right_index_is_smaller() {
        let first = [[10.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let second = [[0.5, 0.0, 0.0]];
        let pairs = cross_pairs(
            &first,
            &second,
            1.0,
            SpatialBackend::BruteForce,
            &ExecutionContext::default(),
        )
        .expect("valid independent coordinate sets");
        assert_eq!(
            pairs,
            vec![CrossPair {
                first: 1,
                second: 0,
                distance_squared: 0.25
            }]
        );
    }

    #[test]
    fn dense_sets_do_not_emit_within_set_pairs() {
        let first = [[0.0; 3]; 4];
        let second = [[0.0; 3]; 3];
        let pairs = cross_pairs(
            &first,
            &second,
            1.0,
            SpatialBackend::CellList,
            &ExecutionContext::default(),
        )
        .expect("valid dense coordinate sets");
        assert_eq!(pairs.len(), first.len() * second.len());
        assert!(pairs.iter().all(|pair| pair.first < 4 && pair.second < 3));
    }
}
