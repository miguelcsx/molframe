//! Implicit, bucketed, cache-friendly three-dimensional k-d tree.
//!
//! # Layout
//!
//! The tree is a perfect binary tree over *index medians*: a node covering the
//! slot range `lo..hi` splits at `lo + (hi - lo) / 2` on the axis of widest
//! extent, after `select_nth_unstable` has placed the median there. Both the
//! shape and every child range therefore follow from the point count alone, so
//! no child pointers are stored. Per node the tree keeps one split value and
//! one axis (heap order, children `2i + 1` and `2i + 2`); per slot it keeps
//! structure-of-arrays coordinates and the atom index. Leaves are buckets of
//! `LEAF_MAX / 2 ..= LEAF_MAX` points scanned linearly, which the compiler
//! vectorises. Splitting by count rather than by value keeps the depth at
//! `ceil(log2(n / LEAF_MAX))` for duplicates, collinear layouts and outliers.
//!
//! Left subtrees hold coordinates `<=` the split value and right subtrees
//! `>=` it, so equal coordinates may sit on either side; searches treat a
//! plane at distance zero as intersecting both.
//!
//! # Concurrency
//!
//! Construction splits the top levels serially and the independent subtrees on
//! the context's workers. Queries are independent and run in blocks on the
//! same workers; no query allocates (scratch is owned per block).

use crate::brute::finite;
use crate::{KdPeriodicOptions, NeighborPair, PeriodicBox, SpatialError};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use std::sync::Mutex;

#[path = "kd_batch.rs"]
mod batch;
#[path = "kd_heap.rs"]
mod heap;
#[path = "kd_periodic.rs"]
mod periodic_helpers;
#[path = "kd_search.rs"]
mod search;
#[path = "kd_stream.rs"]
mod stream;

pub use batch::{KnnOptions, KnnResult, NO_NEIGHBOR};
use periodic_helpers::{for_each_shift, periodic_image_limits, validate_image_budget};

/// Largest leaf bucket. Leaves hold `LEAF_MAX / 2 ..= LEAF_MAX` points.
pub(super) const LEAF_MAX: usize = 32;

#[derive(Clone, Copy, Debug)]
struct Entry {
    position: [f32; 3],
    atom: u32,
}

/// A balanced static k-dimensional tree over selected atoms.
///
/// Construction uses linear-time median partitioning at each level, for
/// O(n log n) total work with `n / LEAF_MAX` split nodes and four flat
/// per-point arrays.
#[derive(Debug)]
pub struct KdTree<'a> {
    positions: &'a [[f32; 3]],
    xs: Vec<f32>,
    ys: Vec<f32>,
    zs: Vec<f32>,
    ids: Vec<u32>,
    split_value: Vec<f32>,
    split_axis: Vec<u8>,
    bounds: Option<([f32; 3], [f32; 3])>,
    periodic: Option<PeriodicBox>,
    periodic_options: KdPeriodicOptions,
    reservation: Option<molframe_core::MemoryReservation>,
}

impl<'a> KdTree<'a> {
    /// Builds a balanced tree over finite target positions.
    ///
    /// # Errors
    ///
    /// Returns an error when a target index is outside `positions`.
    pub fn build(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        periodic: Option<PeriodicBox>,
    ) -> Result<Self, SpatialError> {
        Self::build_with_options(positions, targets, periodic, KdPeriodicOptions::default())
    }

    /// Builds a balanced tree with an explicit periodic image budget.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid options or target indices.
    pub fn build_with_options(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        periodic: Option<PeriodicBox>,
        periodic_options: KdPeriodicOptions,
    ) -> Result<Self, SpatialError> {
        Self::build_parallel(positions, targets, periodic, periodic_options, None)
    }

    pub(super) fn build_parallel(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        periodic: Option<PeriodicBox>,
        periodic_options: KdPeriodicOptions,
        context: Option<&ExecutionContext>,
    ) -> Result<Self, SpatialError> {
        validate_indices(targets, positions.len())?;
        periodic_options.validate()?;

        let mut entries = finite_entries(positions, targets, periodic.as_ref());
        let bounds = bounds_of(&entries);
        let depth = depth_for(entries.len());
        let node_count = (1_usize << depth) - 1;
        let mut split_value = vec![0.0_f32; node_count];
        let mut split_axis = vec![0_u8; node_count];
        for (node, axis, value) in build_splits(&mut entries, depth, context)? {
            split_value[node] = value;
            split_axis[node] = axis;
        }
        let mut xs = Vec::with_capacity(entries.len());
        let mut ys = Vec::with_capacity(entries.len());
        let mut zs = Vec::with_capacity(entries.len());
        let mut ids = Vec::with_capacity(entries.len());
        for entry in &entries {
            xs.push(entry.position[0]);
            ys.push(entry.position[1]);
            zs.push(entry.position[2]);
            ids.push(entry.atom);
        }

        Ok(Self {
            positions,
            xs,
            ys,
            zs,
            ids,
            split_value,
            split_axis,
            bounds,
            periodic,
            periodic_options,
            reservation: None,
        })
    }

    /// Number of finite indexed points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Whether the tree indexes no finite point.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Bytes retained by the finished tree.
    pub(crate) fn retained_bytes(&self) -> usize {
        self.ids.capacity() * (3 * size_of::<f32>() + size_of::<u32>())
            + self.split_value.capacity() * (size_of::<f32>() + size_of::<u8>())
    }

    /// Bytes a build over `targets` may hold at its peak (entries plus result).
    pub(crate) fn build_bytes(targets: usize) -> Option<usize> {
        let depth = depth_for(targets);
        targets
            .checked_mul(size_of::<Entry>() + 3 * size_of::<f32>() + size_of::<u32>())?
            .checked_add(((1_usize << depth) - 1) * (size_of::<f32>() + size_of::<u8>()))
    }

    /// Finds all indexed pairs within `cutoff`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid cutoff or query index.
    pub fn pairs(&self, query: &[u32], cutoff: f32) -> Result<Vec<NeighborPair>, SpatialError> {
        validate_cutoff(cutoff)?;
        validate_indices(query, self.positions.len())?;
        let mut found = Vec::new();
        let mut collector = self.radius_collector(cutoff)?;
        for &atom in query {
            collector.collect(self, atom, |target, squared| {
                found.push(NeighborPair::new(atom, target, squared));
            })?;
        }
        crate::brute::canonicalise(&mut found);
        Ok(found)
    }

    /// Returns the nearest `count` targets in increasing distance order.
    ///
    /// Distances are squared. Equal distances order by atom index.
    ///
    /// # Errors
    ///
    /// Returns an error when `query` is outside the coordinate array.
    pub fn k_nearest(&self, query: u32, count: usize) -> Result<Vec<(u32, f32)>, SpatialError> {
        let index = usize::try_from(query).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let Some(point) = self.positions.get(index).copied() else {
            return Err(SpatialError::AtomOutOfBounds(query));
        };
        if count == 0 || !finite(point) {
            return Ok(Vec::new());
        }
        let mut heap = heap::KnnHeap::new(count);
        self.knn_into(&mut heap, point, Some(query), f32::INFINITY)?;
        Ok(heap
            .sorted()
            .into_iter()
            .map(|(squared, atom)| (atom, squared))
            .collect())
    }
}

/// Returns finite target entries in the tree coordinate system.
///
/// Runtime is `O(T)` with `O(F)` output space for `F` finite targets.
fn finite_entries(
    positions: &[[f32; 3]],
    targets: &[u32],
    periodic: Option<&PeriodicBox>,
) -> Vec<Entry> {
    let mut entries = Vec::with_capacity(targets.len());
    entries.extend(targets.iter().filter_map(|&atom| {
        let index = usize::try_from(atom).ok()?;
        let position = positions.get(index).copied()?;
        if !finite(position) {
            return None;
        }
        let position = match periodic {
            Some(periodic) => periodic.wrap(position),
            None => position,
        };
        Some(Entry { position, atom })
    }));
    entries
}

fn bounds_of(entries: &[Entry]) -> Option<([f32; 3], [f32; 3])> {
    let first = entries.first()?;
    let mut low = first.position;
    let mut high = first.position;
    for entry in entries {
        for axis in 0..3 {
            low[axis] = low[axis].min(entry.position[axis]);
            high[axis] = high[axis].max(entry.position[axis]);
        }
    }
    Some((low, high))
}

/// Smallest depth whose leaves all hold at most `LEAF_MAX` points.
fn depth_for(count: usize) -> usize {
    let mut depth = 0;
    while count.div_ceil(1 << depth) > LEAF_MAX {
        depth += 1;
    }
    depth
}

type Split = (usize, u8, f32);

/// Splits one node's slice at its count median on the widest axis.
fn split_slice(slice: &mut [Entry]) -> (u8, f32, usize) {
    let mut low = [f32::INFINITY; 3];
    let mut high = [f32::NEG_INFINITY; 3];
    for entry in slice.iter() {
        for axis in 0..3 {
            low[axis] = low[axis].min(entry.position[axis]);
            high[axis] = high[axis].max(entry.position[axis]);
        }
    }
    let mut axis = 0_u8;
    for candidate in 1..3_u8 {
        let (current, next) = (usize::from(axis), usize::from(candidate));
        if high[next] - low[next] > high[current] - low[current] {
            axis = candidate;
        }
    }
    let middle = slice.len() / 2;
    let key = usize::from(axis);
    slice.select_nth_unstable_by(middle, |left, right| {
        left.position[key].total_cmp(&right.position[key])
    });
    (axis, slice[middle].position[key], middle)
}

type Task<'s> = Mutex<Option<(usize, &'s mut [Entry])>>;

fn build_subtree(slice: &mut [Entry], node: usize, remaining: usize, out: &mut Vec<Split>) {
    if remaining == 0 {
        return;
    }
    let (axis, value, middle) = split_slice(slice);
    out.push((node, axis, value));
    let (left, right) = slice.split_at_mut(middle);
    build_subtree(left, 2 * node + 1, remaining - 1, out);
    build_subtree(right, 2 * node + 2, remaining - 1, out);
}

/// Computes every split node, parallelising the independent subtrees.
fn build_splits(
    entries: &mut [Entry],
    depth: usize,
    context: Option<&ExecutionContext>,
) -> Result<Vec<Split>, SpatialError> {
    let workers = context.map_or(1, ExecutionContext::worker_budget);
    let parallel = workers > 1 && entries.len() >= 16_384 && depth > 0;
    let mut splits = Vec::new();
    if !parallel {
        build_subtree(entries, 0, depth, &mut splits);
        return Ok(splits);
    }
    // Enough subtrees to balance dynamic scheduling, few enough that the
    // serial top levels stay a small fraction of the work.
    let top = depth.min((workers * 4).next_power_of_two().trailing_zeros() as usize);
    let mut frontier = vec![(0_usize, entries)];
    for _ in 0..top {
        let mut next = Vec::with_capacity(frontier.len() * 2);
        for (node, slice) in frontier {
            let (axis, value, middle) = split_slice(slice);
            splits.push((node, axis, value));
            let (left, right) = slice.split_at_mut(middle);
            next.push((2 * node + 1, left));
            next.push((2 * node + 2, right));
        }
        frontier = next;
    }
    let remaining = depth - top;
    let tasks: Vec<Task> = frontier
        .into_iter()
        .map(|task| Mutex::new(Some(task)))
        .collect();
    let Some(context) = context else {
        return Err(SpatialError::WorkerPanicked);
    };
    let produced = map_blocks_in(BlockPlan::new(tasks.len(), 1), context, |index, _| {
        let mut out = Vec::new();
        let task = tasks.get(index).and_then(|slot| slot.lock().ok()?.take());
        if let Some((node, slice)) = task {
            build_subtree(slice, node, remaining, &mut out);
        }
        out
    })
    .map_err(|_| SpatialError::WorkerPanicked)?;
    for mut block in produced {
        splits.append(&mut block);
    }
    Ok(splits)
}

/// Validates a finite non-negative cutoff.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_cutoff(cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

/// Validates atom indices against a coordinate-array length.
///
/// Runtime is `O(N)` and requires no allocation.
fn validate_indices(atoms: &[u32], position_count: usize) -> Result<(), SpatialError> {
    for &atom in atoms {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if index >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom));
        }
    }

    Ok(())
}

#[cfg(test)]
#[path = "kd_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "kd_adversarial_tests.rs"]
mod adversarial_tests;
