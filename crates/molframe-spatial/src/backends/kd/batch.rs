//! Parallel batch queries over a k-d tree: nearest neighbours and pair lists.

use super::KdTree;
use crate::backends::errors::{block_error, validate_cutoff, validate_indices};
use crate::brute::canonicalise;
use crate::{NeighborPair, SpatialError};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, try_for_each_block_in};

/// Queries handled per scheduled block.
const QUERY_BLOCK: usize = 1024;
/// Index reported for a neighbour slot that no target filled.
pub const NO_NEIGHBOR: u32 = u32::MAX;

/// Options for a batch nearest-neighbour query.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct KnnOptions {
    /// Largest Euclidean distance reported; farther targets leave the slot empty.
    pub max_distance: Option<f32>,
    /// Skip the query atom itself (and, when periodic, its own images).
    pub exclude_self: bool,
}

/// Row-major `(query_count, k)` nearest-neighbour table.
///
/// Rows are ordered by increasing distance, ties by atom index. Slots no
/// target fills hold [`NO_NEIGHBOR`] and an infinite squared distance.
#[derive(Clone, Debug, PartialEq)]
pub struct KnnResult {
    /// Neighbours per query.
    pub k: usize,
    /// Target atom indices, `query_count * k` entries.
    pub indices: Vec<u32>,
    /// Squared distances, `query_count * k` entries.
    pub distances_squared: Vec<f32>,
}

impl KnnOptions {
    /// Options for a table query: `exclude_self` defaults to true when the
    /// points query themselves and to false when separate queries are given.
    #[must_use]
    pub fn for_queries(
        max_distance: Option<f32>,
        exclude_self: Option<bool>,
        separate_queries: bool,
    ) -> Self {
        let exclude_self = match exclude_self {
            Some(choice) => choice,
            None => !separate_queries,
        };
        Self {
            max_distance,
            exclude_self,
        }
    }
}

impl KnnResult {
    /// Consumes the table into its indices and Euclidean distances.
    #[must_use]
    pub fn into_columns(self) -> (Vec<u32>, Vec<f32>) {
        let mut distances = self.distances_squared;
        for value in &mut distances {
            *value = value.sqrt();
        }
        (self.indices, distances)
    }

    /// Number of query rows.
    #[must_use]
    pub fn rows(&self) -> usize {
        match self.indices.len().checked_div(self.k) {
            Some(rows) => rows,
            None => 0,
        }
    }
}

impl KdTree<'_> {
    /// Finds the `k` nearest targets of every query atom in parallel.
    ///
    /// Results are identical at any worker count. Equal distances order by
    /// atom index. Periodic trees report exact minimum-image distances.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid `max_distance`, a query atom outside
    /// the coordinate array, memory admission failure or cancellation.
    pub fn knn(
        &self,
        query: &[u32],
        k: usize,
        options: KnnOptions,
        context: &ExecutionContext,
    ) -> Result<KnnResult, SpatialError> {
        if let Some(distance) = options.max_distance {
            validate_cutoff(distance)?;
        }
        validate_indices(query, self.positions.len())?;
        let slots = query
            .len()
            .checked_mul(k)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let bytes = slots
            .checked_mul(size_of::<u32>() + size_of::<f32>())
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let _reservation = context.try_reserve(bytes)?;
        let mut indices = Vec::new();
        let mut distances = Vec::new();
        indices
            .try_reserve_exact(slots)
            .map_err(|_| SpatialError::NumericRangeExceeded)?;
        distances
            .try_reserve_exact(slots)
            .map_err(|_| SpatialError::NumericRangeExceeded)?;
        if k == 0 || query.is_empty() {
            return Ok(KnnResult {
                k,
                indices,
                distances_squared: distances,
            });
        }
        self.searcher_check(k, options)?;

        let plan = BlockPlan::new(query.len(), QUERY_BLOCK);
        try_for_each_block_in(
            plan,
            context,
            QUERY_BLOCK * k * (size_of::<u32>() + size_of::<f32>()),
            |_, range| {
                let mut searcher = self.knn_searcher(k, options.max_distance)?;
                let chunk = match query.get(range) {
                    Some(chunk) => chunk,
                    None => &[],
                };
                let mut ids = Vec::with_capacity(chunk.len() * k);
                let mut squared = Vec::with_capacity(chunk.len() * k);
                for &atom in chunk {
                    searcher.search(self, atom, options.exclude_self, options.max_distance);
                    let found = searcher.heap.sort_in_place();
                    let limit = options.max_distance.map_or(f32::INFINITY, |d| d * d);
                    let mut filled = 0;
                    for &(distance, id) in found.iter().take(k) {
                        if distance > limit {
                            break;
                        }
                        ids.push(id);
                        squared.push(distance);
                        filled += 1;
                    }
                    for _ in filled..k {
                        ids.push(NO_NEIGHBOR);
                        squared.push(f32::INFINITY);
                    }
                }
                Ok((ids, squared))
            },
            |(ids, squared): (Vec<u32>, Vec<f32>)| {
                indices.extend_from_slice(&ids);
                distances.extend_from_slice(&squared);
                Ok(())
            },
        )
        .map_err(|error| block_error(&error))?;
        Ok(KnnResult {
            k,
            indices,
            distances_squared: distances,
        })
    }

    /// Fails early when the periodic image budget cannot serve the query.
    fn searcher_check(&self, k: usize, options: KnnOptions) -> Result<(), SpatialError> {
        self.knn_searcher(k, options.max_distance).map(|_| ())
    }

    /// Enumerates unique pairs within `cutoff` for the query atoms, in parallel.
    ///
    /// With `same_selection`, the tree must index the same atoms as `query`;
    /// each pair is then produced once, from its smaller atom, and the output
    /// is globally sorted without a global sort when `query` ascends. Cross
    /// selections go through [`canonicalise`] so overlapping selections stay
    /// duplicate-free.
    pub(crate) fn pairs_in(
        &self,
        query: &[u32],
        cutoff: f32,
        same_selection: bool,
        sorted: bool,
        context: &ExecutionContext,
    ) -> Result<Vec<NeighborPair>, SpatialError> {
        validate_cutoff(cutoff)?;
        validate_indices(query, self.positions.len())?;
        // Surface image-budget errors before scheduling any work.
        self.radius_collector(cutoff)?;

        let block = if self.periodic.is_some() { 256 } else { 2048 };
        let plan = BlockPlan::new(query.len(), block);
        let expected = 32 * size_of::<NeighborPair>();
        let mut found: Vec<NeighborPair> = Vec::new();
        try_for_each_block_in(
            plan,
            context,
            block.saturating_mul(expected),
            |_, range| {
                let mut collector = self.radius_collector(cutoff)?;
                let mut out = Vec::new();
                let mut hits: Vec<(u32, f32)> = Vec::new();
                let chunk = match query.get(range) {
                    Some(chunk) => chunk,
                    None => &[],
                };
                for &atom in chunk {
                    hits.clear();
                    collector.collect(self, atom, |target, squared| {
                        if !same_selection || target > atom {
                            hits.push((target, squared));
                        }
                    })?;
                    if sorted || same_selection {
                        hits.sort_unstable_by_key(|hit| hit.0);
                    }
                    out.extend(
                        hits.iter()
                            .map(|&(target, squared)| NeighborPair::new(atom, target, squared)),
                    );
                }
                Ok(out)
            },
            |mut block: Vec<NeighborPair>| {
                found.append(&mut block);
                Ok(())
            },
        )
        .map_err(|error| block_error(&error))?;

        let ascending = query.windows(2).all(|pair| pair[0] < pair[1]);
        if !same_selection || (sorted && !ascending) {
            canonicalise(&mut found);
        }
        Ok(found)
    }
}
