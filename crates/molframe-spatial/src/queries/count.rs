//! Scalar contact counts with bounded, ordered partial reductions.

use super::{PairQuery, indices::IndexWorkspace, search::validate_cutoff, visit};
use crate::{CellList, SpatialBackend, SpatialError};
use molframe_core::parallel::{BlockExecutionError, BlockPlan, try_for_each_block_in};

/// Counts unique unordered pairs without storing pairs or all block results.
///
/// Cell queries share one charged index across dynamically scheduled blocks.
/// Fixed 64-cell or 256-query boundaries are independent of worker count. Other
/// backends use their bounded visitor. Memory scales with the index and the
/// admitted worker window, independently of the number of matching pairs.
///
/// # Errors
///
/// Returns input, memory, cancellation or worker errors before returning a count.
pub fn count_pairs_within(query: &PairQuery<'_>) -> Result<u64, SpatialError> {
    validate_cutoff(query.cutoff)?;
    let indices = IndexWorkspace::new(
        query.left,
        query.right,
        query.positions.len(),
        query.context,
    )?;
    let left = indices.left();
    let right = indices.right();
    let plan = query.options.plan(left.len(), right.len(), query.cutoff)?;
    let mut total = 0_u64;
    if !matches!(
        plan.backend,
        SpatialBackend::CellList | SpatialBackend::NeighborList
    ) {
        visit::visit_indices(query, left, right, |_| total += 1)?;
        return Ok(total);
    }
    let index = CellList::build_in(
        query.positions,
        right,
        query.cutoff,
        query.periodic.copied(),
        query.options.cell_grid,
        query.context,
    )?;
    let same = left == right;
    let owned_cells = if same { index.owned_cell_count() } else { None };
    let blocks = owned_cells.map_or_else(
        || BlockPlan::new(left.len(), if same { 64 } else { 256 }),
        |count| BlockPlan::new(count, 64),
    );
    try_for_each_block_in(
        blocks,
        query.context,
        size_of::<u64>(),
        |_, range| {
            let mut count = 0_u64;
            if owned_cells.is_some() {
                index.for_each_owned_cell(range, query.cutoff, |_| count += 1)?;
                return Ok(count);
            }
            if same {
                index.for_each_upper_candidate(&left[range], query.cutoff, |_| count += 1)?;
                return Ok(count);
            }
            index.for_each_candidate(&left[range], query.cutoff, |a, b, distance| {
                visit::emit_unique_cross_candidate(a, b, distance, left, right, &mut |_| {
                    count += 1;
                });
            })?;
            Ok(count)
        },
        |count| {
            total = total
                .checked_add(count)
                .ok_or(SpatialError::NumericRangeExceeded)?;
            Ok(())
        },
    )
    .map_err(|error| match error {
        BlockExecutionError::Memory(error) => SpatialError::Memory(error),
        BlockExecutionError::Cancelled => SpatialError::Cancelled,
        BlockExecutionError::Worker(_) => SpatialError::WorkerPanicked,
        BlockExecutionError::Operation(error) => error,
    })?;
    Ok(total)
}

#[cfg(test)]
#[path = "count_tests.rs"]
mod tests;
