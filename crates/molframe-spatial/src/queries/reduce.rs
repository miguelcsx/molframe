//! Fixed-radius queries that reduce pairs instead of returning them.
//!
//! A materialising pair query retains `O(pairs)` bytes. At a billion atoms and
//! a typical contact cutoff that is hundreds of gigabytes for a result the
//! caller immediately collapses into a per-residue map or a union-find. This
//! module folds each pair as it is produced, so peak retained bytes track the
//! accumulator rather than the pair count (FR-516) while keeping the worker
//! parallelism a materialising query already had.

use super::search::validate_cutoff;
use super::visit::{emit_unique_cross_candidate, visit_indices};
use crate::{
    CellList, NeighborPair, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions,
};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use molframe_core::selection::AtomSelection;

/// Everything a fixed-radius query needs, so a reduction takes three arguments.
#[derive(Clone, Copy, Debug)]
pub struct PairQuery<'a> {
    /// Coordinates that every index in both selections refers to.
    pub positions: &'a [[f32; 3]],
    /// The selection each pair's first endpoint is drawn from.
    pub left: &'a AtomSelection,
    /// The selection each pair's second endpoint is drawn from.
    pub right: &'a AtomSelection,
    /// The inclusive distance within which two atoms form a pair.
    pub cutoff: f32,
    /// Backend and grid selection.
    pub options: SpatialSearchOptions,
    /// The periodic cell to apply the minimum image convention over, if any.
    pub periodic: Option<&'a PeriodicBox>,
    /// Shared execution budget and native worker pool.
    pub context: &'a ExecutionContext,
}

/// Folds every fixed-radius pair into per-block accumulators.
///
/// The accumulators come back in block order, which is a function of the cell
/// count and never of `workers`, so merging them in the order returned gives
/// the same result at any worker count (FR-515). The caller supplies the merge,
/// because only it knows whether the reduction is a minimum, a union or a sum.
///
/// Falls back to a single accumulator when the workload has no block
/// decomposition — a cross-selection query, a non-cell backend, a periodic box
/// or one worker — so the caller's merge is written once either way.
///
/// # Errors
///
/// Returns the same validation and backend errors as the materialising query.
pub fn reduce_pairs_within_unsorted<T, I, F>(
    query: &PairQuery<'_>,
    init: I,
    fold: F,
) -> Result<Vec<T>, SpatialError>
where
    T: Send,
    I: Fn() -> T + Sync,
    F: Fn(&mut T, NeighborPair) + Sync,
{
    let &PairQuery {
        positions,
        left,
        right,
        cutoff,
        options,
        periodic,
        context,
    } = query;

    validate_cutoff(cutoff)?;
    let indices = super::indices::IndexWorkspace::new(left, right, positions.len(), context)?;
    let left_indices = indices.left();
    let right_indices = indices.right();
    let plan = options.plan(
        left_indices.len(),
        right_indices.len(),
        periodic.is_some(),
        cutoff,
    )?;

    // A cell list over a non-periodic workload is the only shape with a block
    // decomposition: cells own disjoint same-selection pairs, and each query
    // atom owns its own cross-selection pairs.
    let workers = context.worker_budget();
    let blockable = plan.backend == SpatialBackend::CellList && periodic.is_none() && workers > 1;
    let same_selection = left_indices == right_indices;

    if blockable {
        let index = CellList::build_in(
            positions,
            right_indices,
            cutoff,
            periodic,
            options.cell_grid,
            context,
        )?;

        if same_selection {
            if let Some(parts) = crate::backends::cell::reduce_pairs_same_selection_parallel(
                &index,
                left_indices,
                cutoff,
                context,
                &init,
                &fold,
            )? {
                return Ok(parts);
            }
        } else {
            return reduce_cross_selection(
                &index,
                left_indices,
                right_indices,
                cutoff,
                context,
                &init,
                &fold,
            );
        }
    }

    let mut single = init();
    visit_indices(query, left_indices, right_indices, |pair| {
        fold(&mut single, pair);
    })?;
    Ok(vec![single])
}

#[cfg(test)]
#[path = "reduce_tests.rs"]
mod tests;

/// Folds cross-selection pairs into per-block accumulators.
///
/// Every pair belongs to exactly one query atom, so partitioning the query
/// indices partitions the pairs. Block boundaries come from the query count and
/// never from `workers`, so merging in the order returned is worker-count
/// independent (FR-515).
fn reduce_cross_selection<T, I, F>(
    index: &CellList<'_>,
    left: &[u32],
    right: &[u32],
    cutoff: f32,
    context: &ExecutionContext,
    init: &I,
    fold: &F,
) -> Result<Vec<T>, SpatialError>
where
    T: Send,
    I: Fn() -> T + Sync,
    F: Fn(&mut T, NeighborPair) + Sync,
{
    /// Query atoms per block. Local density varies, so blocks stay well below
    /// one worker's share rather than being one share each.
    const BLOCK_QUERIES: usize = 256;

    let plan = BlockPlan::new(left.len(), BLOCK_QUERIES);
    let produced = map_blocks_in(plan, context, |_, range| {
        let mut accumulator = init();
        let Some(chunk) = left.get(range) else {
            return Ok(accumulator);
        };
        index.for_each_candidate(chunk, cutoff, |left_atom, right_atom, squared| {
            emit_unique_cross_candidate(left_atom, right_atom, squared, left, right, &mut |pair| {
                fold(&mut accumulator, pair);
            });
        })?;
        Ok(accumulator)
    })
    .map_err(|_| SpatialError::WorkerPanicked)?;

    produced.into_iter().collect()
}
