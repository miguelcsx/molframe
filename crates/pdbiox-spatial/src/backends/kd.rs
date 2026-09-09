//! Balanced three-dimensional k-d tree.

use crate::brute::{canonicalise, distance_squared, finite};
use crate::{KdPeriodicOptions, NeighborPair, PeriodicBox, SpatialError};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

#[path = "kd_periodic.rs"]
mod periodic_helpers;
#[path = "kd_stream.rs"]
mod stream;
use periodic_helpers::{for_each_shift, periodic_image_limits, validate_image_budget};

#[derive(Clone, Copy, Debug)]
struct Node {
    atom: u32,
    position: [f32; 3],
    axis: usize,
    left: Option<usize>,
    right: Option<usize>,
}

/// A balanced static k-dimensional tree over selected atoms.
///
/// Construction uses linear-time median partitioning at each level, for
/// O(n log n) total work and one flat node allocation.
#[derive(Debug)]
pub struct KdTree<'a> {
    positions: &'a [[f32; 3]],
    nodes: Vec<Node>,
    root: Option<usize>,
    periodic: Option<&'a PeriodicBox>,
    periodic_options: KdPeriodicOptions,
    reservation: Option<pdbiox_core::MemoryReservation>,
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
        periodic: Option<&'a PeriodicBox>,
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
        periodic: Option<&'a PeriodicBox>,
        periodic_options: KdPeriodicOptions,
    ) -> Result<Self, SpatialError> {
        validate_indices(targets, positions.len())?;
        periodic_options.validate()?;

        let mut entries = finite_entries(positions, targets, periodic);
        let mut nodes = Vec::with_capacity(entries.len());
        let root = build_nodes(&mut entries, 0, &mut nodes);

        Ok(Self {
            positions,
            nodes,
            root,
            periodic,
            periodic_options,
            reservation: None,
        })
    }

    /// Finds all indexed pairs within `cutoff`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid cutoff or query index.
    pub fn pairs(&self, query: &[u32], cutoff: f32) -> Result<Vec<NeighborPair>, SpatialError> {
        validate_cutoff(cutoff)?;
        validate_indices(query, self.positions.len())?;

        let cutoff_squared = cutoff * cutoff;
        let mut found = Vec::new();

        for &atom in query {
            let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
            let Some(position) = self.positions.get(index).copied() else {
                continue;
            };

            if !finite(position) {
                continue;
            }

            self.append_radius_pairs(atom, position, cutoff, cutoff_squared, &mut found)?;
        }

        canonicalise(&mut found);
        Ok(found)
    }

    /// Returns the nearest `count` targets in increasing distance order.
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

        let Some(root) = self.root else {
            return Ok(Vec::new());
        };

        if let Some(periodic) = self.periodic {
            return self.periodic_nearest(root, query, point, count, periodic);
        }

        let capacity = count
            .checked_add(1)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let mut heap = BinaryHeap::with_capacity(capacity);

        self.nearest_search(root, query, point, count, &mut heap);

        Ok(sorted_candidates(heap))
    }

    fn append_radius_pairs(
        &self,
        atom: u32,
        position: [f32; 3],
        cutoff: f32,
        cutoff_squared: f32,
        found: &mut Vec<NeighborPair>,
    ) -> Result<(), SpatialError> {
        let Some(root) = self.root else {
            return Ok(());
        };
        let Some(periodic) = self.periodic else {
            self.radius_search(root, position, cutoff_squared, &mut |_, target, squared| {
                if atom != target {
                    found.push(NeighborPair::new(atom, target, squared));
                }
            });
            return Ok(());
        };

        let limits = periodic_image_limits(periodic, cutoff)?;
        validate_image_budget(limits, self.periodic_options)?;
        let wrapped = periodic.wrap(position);
        for_each_shift(limits, |shift| {
            let image = periodic.translated(wrapped, shift)?;
            self.radius_search(root, image, cutoff_squared, &mut |_, target, _| {
                if atom == target {
                    return;
                }
                let Ok(index) = usize::try_from(target) else {
                    return;
                };
                let Some(target_position) = self.positions.get(index).copied() else {
                    return;
                };
                let squared = periodic.distance_squared(position, target_position);
                if squared <= cutoff_squared {
                    found.push(NeighborPair::new(atom, target, squared));
                }
            });
            Ok(())
        })?;
        Ok(())
    }

    /// Recursively visits k-d nodes that can intersect a query sphere.
    ///
    /// The tree is structurally balanced, so recursion depth is `O(log N)`.
    /// No per-query traversal vector is allocated.
    fn radius_search(
        &self,
        node_index: usize,
        point: [f32; 3],
        cutoff_squared: f32,
        visit: &mut impl FnMut(usize, u32, f32),
    ) {
        let Some(node) = self.nodes.get(node_index).copied() else {
            return;
        };

        let position = node.position;

        let squared = distance_squared(point, position, None);

        if squared <= cutoff_squared {
            visit(node_index, node.atom, squared);
        }

        let delta = point[node.axis] - position[node.axis];
        let (near, far) = ordered_children(node, delta);

        if let Some(near) = near {
            self.radius_search(near, point, cutoff_squared, visit);
        }

        if delta * delta <= cutoff_squared
            && let Some(far) = far
        {
            self.radius_search(far, point, cutoff_squared, visit);
        }
    }

    /// Recursively performs branch-and-bound nearest-neighbour search.
    ///
    /// The near child is processed before the far-child pruning decision,
    /// tightening the heap radius as early as possible.
    fn nearest_search(
        &self,
        node_index: usize,
        query: u32,
        point: [f32; 3],
        count: usize,
        heap: &mut BinaryHeap<Candidate>,
    ) {
        let Some(node) = self.nodes.get(node_index).copied() else {
            return;
        };

        let position = node.position;

        let squared = distance_squared(point, position, None);

        if node.atom != query {
            push_candidate(heap, Candidate::new(node.atom, squared), count);
        }

        let axis = node.axis;
        let delta = point[axis] - position[axis];
        let (near, far) = ordered_children(node, delta);

        if let Some(near) = near {
            self.nearest_search(near, query, point, count, heap);
        }

        if should_visit_far(heap, count, delta)
            && let Some(far) = far
        {
            self.nearest_search(far, query, point, count, heap);
        }
    }

    fn periodic_nearest(
        &self,
        root: usize,
        query: u32,
        point: [f32; 3],
        count: usize,
        periodic: &PeriodicBox,
    ) -> Result<Vec<(u32, f32)>, SpatialError> {
        const NEAREST_IMAGE_LIMITS: [i64; 3] = [2; 3];
        validate_image_budget(NEAREST_IMAGE_LIMITS, self.periodic_options)?;
        let wrapped = periodic.wrap(point);
        let mut best_by_atom = BTreeMap::new();

        for_each_shift(NEAREST_IMAGE_LIMITS, |shift| {
            let image = periodic.translated(wrapped, shift)?;
            let capacity = count
                .checked_add(1)
                .ok_or(SpatialError::NumericRangeExceeded)?;
            let mut image_heap = BinaryHeap::with_capacity(capacity);
            self.nearest_search(root, query, image, count, &mut image_heap);
            for (atom, _) in sorted_candidates(image_heap) {
                let index =
                    usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
                let Some(target) = self.positions.get(index).copied() else {
                    continue;
                };
                let squared = periodic.distance_squared(point, target);
                best_by_atom
                    .entry(atom)
                    .and_modify(|best: &mut f32| *best = best.min(squared))
                    .or_insert(squared);
            }
            Ok(())
        })?;

        let capacity = count
            .checked_add(1)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let mut heap = BinaryHeap::with_capacity(capacity);
        for (atom, squared) in best_by_atom {
            push_candidate(&mut heap, Candidate::new(atom, squared), count);
        }
        Ok(sorted_candidates(heap))
    }
}

#[derive(Clone, Copy, Debug)]
struct Entry {
    atom: u32,
    position: [f32; 3],
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
        Some(Entry { atom, position })
    }));
    entries
}

/// Builds one balanced k-d subtree by median partitioning.
///
/// Construction satisfies `O(N log N)` total expected/average partition work
/// and `O(log N)` recursion depth because each call splits by index median.
fn build_nodes(entries: &mut [Entry], depth: usize, nodes: &mut Vec<Node>) -> Option<usize> {
    if entries.is_empty() {
        return None;
    }

    let axis = depth % 3;
    let middle = entries.len() / 2;

    entries.select_nth_unstable_by(middle, |left, right| {
        left.position[axis]
            .total_cmp(&right.position[axis])
            .then_with(|| left.atom.cmp(&right.atom))
    });

    let (left, rest) = entries.split_at_mut(middle);
    let (median, right) = rest.split_first_mut()?;
    let entry = *median;

    let node_index = nodes.len();
    nodes.push(Node {
        atom: entry.atom,
        position: entry.position,
        axis,
        left: None,
        right: None,
    });

    let left_node = build_nodes(left, depth + 1, nodes);
    let right_node = build_nodes(right, depth + 1, nodes);

    if let Some(node) = nodes.get_mut(node_index) {
        node.left = left_node;
        node.right = right_node;
    }

    Some(node_index)
}

/// Orders a node's children by the query's position relative to its split.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn ordered_children(node: Node, delta: f32) -> (Option<usize>, Option<usize>) {
    if delta <= 0.0 {
        (node.left, node.right)
    } else {
        (node.right, node.left)
    }
}

/// Returns whether branch-and-bound search must inspect the far subtree.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn should_visit_far(heap: &BinaryHeap<Candidate>, count: usize, delta: f32) -> bool {
    if heap.len() < count {
        return true;
    }

    heap.peek()
        .is_some_and(|candidate| delta * delta <= candidate.squared)
}

#[derive(Clone, Copy, Debug)]
struct Candidate {
    atom: u32,
    squared: f32,
}

impl Candidate {
    /// Creates a nearest-neighbour heap candidate.
    ///
    /// Runtime and auxiliary space are `O(1)`.
    const fn new(atom: u32, squared: f32) -> Self {
        Self { atom, squared }
    }
}

impl PartialEq for Candidate {
    /// Tests exact candidate identity using the floating-point bit pattern.
    fn eq(&self, other: &Self) -> bool {
        self.atom == other.atom && self.squared.to_bits() == other.squared.to_bits()
    }
}

impl Eq for Candidate {}

impl PartialOrd for Candidate {
    /// Delegates candidate ordering to its total ordering.
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Candidate {
    /// Orders the heap by distance and then atom index.
    fn cmp(&self, other: &Self) -> Ordering {
        self.squared
            .total_cmp(&other.squared)
            .then_with(|| self.atom.cmp(&other.atom))
    }
}

/// Inserts a candidate into a fixed-size max heap.
///
/// Runtime is `O(log K)` and heap space remains `O(K)`.
fn push_candidate(heap: &mut BinaryHeap<Candidate>, candidate: Candidate, count: usize) {
    heap.push(candidate);

    if heap.len() > count {
        let _discarded = heap.pop();
    }
}

/// Converts the bounded heap to increasing distance order.
///
/// Runtime is `O(K log K)` as provided by `BinaryHeap::into_sorted_vec`.
fn sorted_candidates(heap: BinaryHeap<Candidate>) -> Vec<(u32, f32)> {
    heap.into_sorted_vec()
        .into_iter()
        .map(|candidate| (candidate.atom, candidate.squared))
        .collect()
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
