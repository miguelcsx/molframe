//! TM-scores for many structure pairs, parallel across pairs.
//!
//! A single [`crate::tm_score`] already spreads its seeds over the shared
//! pool, which is the wrong grain for thousands of small pairs: each call
//! pays fork-join overhead and allocates its own scratch. Here whole pairs are
//! the unit of parallel work and each pair runs its seeds inline, so the
//! values equal the single-pair function bit for bit at any worker count.

use crate::failure::CompareError;
use crate::seed_search::Schedule;
use crate::superposed::tm_score_scheduled;
use crate::tm_screen::{Stage, Tally};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockExecutionError, BlockPlan, try_for_each_block_in};
use molframe_spatial::SpatialError;

/// Pairs scored per scheduled block.
const PAIRS_PER_BLOCK: usize = 4;

/// Scratch bytes charged per in-flight pair and coordinate.
const BYTES_PER_COORDINATE: usize = 96;

/// Scores `model = structures[a]` against `reference = structures[b]` for each
/// `(a, b)` in `pairs`, in the order given.
///
/// Every pair needs equal lengths, as for [`crate::tm_score`]. Memory is
/// `8 * pairs.len()` for the result plus scratch for a bounded window of pairs.
///
/// # Errors
///
/// Returns [`CompareError::InvalidMapping`] for an out-of-range structure index,
/// the first scoring error in pair order, or a memory, cancellation or worker
/// failure.
pub fn tm_score_pairs(
    structures: &[&[[f32; 3]]],
    pairs: &[(usize, usize)],
    context: &ExecutionContext,
) -> Result<Vec<f64>, CompareError> {
    if pairs
        .iter()
        .any(|&(a, b)| a >= structures.len() || b >= structures.len())
    {
        return Err(CompareError::InvalidMapping);
    }
    run(
        structures,
        pairs.len(),
        context,
        |index| pairs.get(index).copied(),
        0,
        |&(a, b), (): &mut ()| {
            tm_score_scheduled(structures[a], structures[b], Schedule::Inline)
                .map(|score| (score, Stage::Scored))
        },
    )
    .map(|(scores, _)| scores)
}

/// Scores every unordered pair of `structures`, returning the condensed upper
/// triangle `(0,1), (0,2), ..., (n-2,n-1)` with `structures[i]` as the model
/// and `structures[j]` as the reference for `i < j`.
///
/// The result holds `n (n - 1) / 2` values, so ten thousand structures need
/// 400 MB; the cost is that many full seed searches, which is why this is the
/// exhaustive reference to filter against rather than a screening step.
///
/// # Errors
///
/// As for [`tm_score_pairs`], without the index error.
pub fn tm_score_all_pairs(
    structures: &[&[[f32; 3]]],
    context: &ExecutionContext,
) -> Result<Vec<f64>, CompareError> {
    let count = structures.len();
    let total = count.checked_mul(count.saturating_sub(1)).map(|v| v / 2);
    let total = total.ok_or(CompareError::InvalidMapping)?;
    run(
        structures,
        total,
        context,
        |index| condensed_pair(count, index),
        0,
        |&(a, b), (): &mut ()| {
            tm_score_scheduled(structures[a], structures[b], Schedule::Inline)
                .map(|score| (score, Stage::Scored))
        },
    )
    .map(|(scores, _)| scores)
}

/// Maps a condensed upper-triangle index back to `(row, column)`.
pub(crate) fn condensed_pair(count: usize, index: usize) -> Option<(usize, usize)> {
    let before = |row: usize| row * (2 * count - row - 1) / 2;
    // The row is the last one whose first entry does not exceed `index`.
    let (mut low, mut high) = (0, count.saturating_sub(1));
    while low < high {
        let middle = low + (high - low).div_ceil(2);
        if before(middle) <= index {
            low = middle;
        } else {
            high = middle - 1;
        }
    }
    (low + 1 < count && index < before(low + 1)).then(|| (low, low + 1 + index - before(low)))
}

/// Scores `total` pairs in blocks, giving `score` a workspace per block.
///
/// `extra_bytes` is the workspace size a block charges beyond the coordinates.
pub(crate) fn run<W: Default>(
    structures: &[&[[f32; 3]]],
    total: usize,
    context: &ExecutionContext,
    pair_at: impl Fn(usize) -> Option<(usize, usize)> + Sync,
    extra_bytes: usize,
    score: impl Fn(&(usize, usize), &mut W) -> Result<(f64, Stage), CompareError> + Sync,
) -> Result<(Vec<f64>, Tally), CompareError> {
    let mut scores = Vec::new();
    scores
        .try_reserve_exact(total)
        .map_err(|_| CompareError::Spatial(SpatialError::NumericRangeExceeded))?;
    let longest = match structures.iter().map(|s| s.len()).max() {
        Some(longest) => longest,
        None => 0,
    };
    let bytes_per_block = PAIRS_PER_BLOCK
        .saturating_mul(longest)
        .saturating_mul(BYTES_PER_COORDINATE)
        .saturating_add(extra_bytes);
    let mut tally = Tally::default();
    try_for_each_block_in(
        BlockPlan::new(total, PAIRS_PER_BLOCK),
        context,
        bytes_per_block,
        |_, range| {
            let mut workspace = W::default();
            let mut block = Tally::default();
            let mut values = Vec::with_capacity(range.len());
            for index in range {
                let pair = pair_at(index).ok_or(CompareError::InvalidMapping)?;
                let (value, stage) = score(&pair, &mut workspace)?;
                block.record(stage);
                values.push(value);
            }
            Ok((values, block))
        },
        |(mut values, block)| {
            scores.append(&mut values);
            tally.add(&block);
            Ok(())
        },
    )
    .map_err(|error| match error {
        BlockExecutionError::Memory(error) => CompareError::Spatial(SpatialError::Memory(error)),
        BlockExecutionError::Cancelled => CompareError::Spatial(SpatialError::Cancelled),
        BlockExecutionError::Worker(_) => CompareError::Spatial(SpatialError::WorkerPanicked),
        BlockExecutionError::Operation(error) => error,
    })?;
    Ok((scores, tally))
}

#[cfg(test)]
#[path = "tm_batch_tests.rs"]
mod tests;
