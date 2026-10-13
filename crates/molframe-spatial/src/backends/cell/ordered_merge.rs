//! Parallel assembly of sorted pair blocks into one canonical pair list.
//!
//! Each block arrives already sorted. The first-atom axis is cut into fixed
//! ranges, so every range owns a disjoint, contiguous slice of the final
//! order. A range gathers its slice from every block by binary search and
//! sorts it independently; ranges are appended in ascending order. The result
//! depends only on the pair set, never on the worker count.

use crate::brute::canonicalise;
use crate::{NeighborPair, SpatialError};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockExecutionError, BlockPlan, try_for_each_block_in};
use std::mem::size_of;

/// First-atom indices covered by one merge range.
const ATOMS_PER_RANGE: usize = 8_192;

/// Maps a block-execution failure onto the spatial error vocabulary.
pub(super) fn spatial_error(error: &BlockExecutionError<SpatialError>) -> SpatialError {
    match error {
        BlockExecutionError::Memory(error) => SpatialError::Memory(*error),
        BlockExecutionError::Cancelled => SpatialError::Cancelled,
        BlockExecutionError::Worker(_) => SpatialError::WorkerPanicked,
        BlockExecutionError::Operation(error) => *error,
    }
}

/// Merges individually sorted blocks into one sorted, duplicate-free vector.
///
/// `atoms` bounds every `first` index. Runtime is `O(P log P / W)` for `P`
/// pairs on `W` workers; auxiliary space is one extra copy of the pairs plus a
/// bounded window of ranges.
pub(super) fn merge_sorted_blocks(
    blocks: &[Vec<NeighborPair>],
    atoms: usize,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    let total: usize = blocks.iter().map(Vec::len).sum();
    let mut merged = Vec::new();
    merged
        .try_reserve_exact(total)
        .map_err(|_| SpatialError::NumericRangeExceeded)?;
    if total == 0 {
        return Ok(merged);
    }

    let ranges = atoms.div_ceil(ATOMS_PER_RANGE).max(1);
    let bytes_per_range = total
        .div_ceil(ranges)
        .saturating_mul(2)
        .saturating_mul(size_of::<NeighborPair>());
    try_for_each_block_in(
        BlockPlan::new(ranges, 1),
        context,
        bytes_per_range,
        |range, _| Ok(gather_range(blocks, range, ranges, atoms)),
        |mut part| {
            merged.append(&mut part);
            Ok(())
        },
    )
    .map_err(|error| spatial_error(&error))?;
    Ok(merged)
}

/// Collects and sorts the pairs whose first atom lies in one range.
fn gather_range(
    blocks: &[Vec<NeighborPair>],
    range: usize,
    ranges: usize,
    atoms: usize,
) -> Vec<NeighborPair> {
    let low = range * ATOMS_PER_RANGE;
    let high = if range + 1 == ranges {
        atoms.max(low)
    } else {
        low + ATOMS_PER_RANGE
    };
    let below = |bound: usize, block: &[NeighborPair]| {
        block.partition_point(|pair| (pair.first as usize) < bound)
    };
    let slices: Vec<&[NeighborPair]> = blocks
        .iter()
        .filter_map(|block| block.get(below(low, block)..below(high, block)))
        .filter(|slice| !slice.is_empty())
        .collect();
    let mut part = Vec::with_capacity(slices.iter().map(|slice| slice.len()).sum());
    for slice in slices {
        part.extend_from_slice(slice);
    }
    canonicalise(&mut part);
    part
}
