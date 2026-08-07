//! Balanced three-dimensional k-d tree.

use crate::brute::{canonicalise, distance_squared, finite};
use crate::{NeighborPair, PeriodicBox, SpatialError};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone, Copy, Debug)]
struct Node {
    atom: u32,
    axis: u8,
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
    targets: Vec<u32>,
    nodes: Vec<Node>,
    root: Option<usize>,
    periodic: Option<&'a PeriodicBox>,
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
        for atom in targets {
            if *atom as usize >= positions.len() {
                return Err(SpatialError::AtomOutOfBounds(*atom));
            }
        }
        let mut indices: Vec<u32> = targets
            .iter()
            .copied()
            .filter(|atom| positions.get(*atom as usize).copied().is_some_and(finite))
            .collect();
        let mut nodes = Vec::with_capacity(indices.len());
        let root = if periodic.is_some() {
            None
        } else {
            build_nodes(positions, &mut indices, 0, &mut nodes)
        };
        Ok(Self {
            positions,
            targets: targets.to_vec(),
            nodes,
            root,
            periodic,
        })
    }

    /// Finds all indexed pairs within `cutoff`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid cutoff or query index.
    pub fn pairs(&self, query: &[u32], cutoff: f32) -> Result<Vec<NeighborPair>, SpatialError> {
        if !cutoff.is_finite() || cutoff < 0.0 {
            return Err(SpatialError::InvalidCutoff);
        }
        for atom in query {
            if *atom as usize >= self.positions.len() {
                return Err(SpatialError::AtomOutOfBounds(*atom));
            }
        }
        if self.periodic.is_some() {
            return Ok(crate::brute::pairs(
                self.positions,
                query,
                &self.targets,
                cutoff * cutoff,
                self.periodic,
            ));
        }
        let cutoff_squared = cutoff * cutoff;
        let mut found = Vec::new();
        for atom in query {
            let Some(position) = self.positions.get(*atom as usize).copied() else {
                continue;
            };
            if finite(position) {
                self.radius(position, cutoff_squared, |target, squared| {
                    if *atom != target {
                        found.push(NeighborPair::new(*atom, target, squared));
                    }
                });
            }
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
        let Some(position) = self.positions.get(query as usize).copied() else {
            return Err(SpatialError::AtomOutOfBounds(query));
        };
        if count == 0 || !finite(position) {
            return Ok(Vec::new());
        }
        if self.periodic.is_some() {
            return Ok(nearest_brute(
                self.positions,
                &self.targets,
                query,
                position,
                count,
                self.periodic,
            ));
        }

        let mut heap = BinaryHeap::with_capacity(count.saturating_add(1));
        let mut stack = Vec::new();
        if let Some(root) = self.root {
            stack.push(root);
        }
        while let Some(node_index) = stack.pop() {
            let Some(node) = self.nodes.get(node_index).copied() else {
                continue;
            };
            let Some(target) = self.positions.get(node.atom as usize).copied() else {
                continue;
            };
            let squared = distance_squared(position, target, None);
            if node.atom != query {
                push_candidate(&mut heap, Candidate::new(node.atom, squared), count);
            }
            let axis = node.axis as usize;
            let delta = position[axis] - target[axis];
            let (near, far) = if delta <= 0.0 {
                (node.left, node.right)
            } else {
                (node.right, node.left)
            };
            if let Some(near) = near {
                stack.push(near);
            }
            let limit = heap
                .peek()
                .map_or(f32::INFINITY, |candidate| candidate.squared);
            if (heap.len() < count || delta * delta <= limit)
                && let Some(far) = far
            {
                stack.push(far);
            }
        }
        Ok(sorted_candidates(heap))
    }

    fn radius(&self, point: [f32; 3], cutoff_squared: f32, mut visit: impl FnMut(u32, f32)) {
        let mut stack = Vec::new();
        if let Some(root) = self.root {
            stack.push(root);
        }
        while let Some(node_index) = stack.pop() {
            let Some(node) = self.nodes.get(node_index).copied() else {
                continue;
            };
            let Some(position) = self.positions.get(node.atom as usize).copied() else {
                continue;
            };
            let squared = distance_squared(point, position, None);
            if squared <= cutoff_squared {
                visit(node.atom, squared);
            }
            let axis = node.axis as usize;
            let delta = point[axis] - position[axis];
            if delta <= 0.0 {
                if let Some(left) = node.left {
                    stack.push(left);
                }
                if delta * delta <= cutoff_squared
                    && let Some(right) = node.right
                {
                    stack.push(right);
                }
            } else {
                if let Some(right) = node.right {
                    stack.push(right);
                }
                if delta * delta <= cutoff_squared
                    && let Some(left) = node.left
                {
                    stack.push(left);
                }
            }
        }
    }
}

fn build_nodes(
    positions: &[[f32; 3]],
    indices: &mut [u32],
    depth: usize,
    nodes: &mut Vec<Node>,
) -> Option<usize> {
    if indices.is_empty() {
        return None;
    }
    let axis = depth % 3;
    let middle = indices.len() / 2;
    indices.select_nth_unstable_by(middle, |left, right| {
        coordinate(positions, *left, axis)
            .total_cmp(&coordinate(positions, *right, axis))
            .then_with(|| left.cmp(right))
    });
    let (left, rest) = indices.split_at_mut(middle);
    let (median, right) = rest.split_first_mut()?;
    let node_index = nodes.len();
    nodes.push(Node {
        atom: *median,
        axis: axis as u8,
        left: None,
        right: None,
    });
    let left_node = build_nodes(positions, left, depth + 1, nodes);
    let right_node = build_nodes(positions, right, depth + 1, nodes);
    if let Some(node) = nodes.get_mut(node_index) {
        node.left = left_node;
        node.right = right_node;
    }
    Some(node_index)
}

fn coordinate(positions: &[[f32; 3]], atom: u32, axis: usize) -> f32 {
    positions
        .get(atom as usize)
        .and_then(|position| position.get(axis))
        .copied()
        .map_or(0.0, |value| value)
}

#[derive(Clone, Copy, Debug)]
struct Candidate {
    atom: u32,
    squared: f32,
}

impl Candidate {
    const fn new(atom: u32, squared: f32) -> Self {
        Self { atom, squared }
    }
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.atom == other.atom && self.squared.to_bits() == other.squared.to_bits()
    }
}

impl Eq for Candidate {}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.squared
            .total_cmp(&other.squared)
            .then_with(|| self.atom.cmp(&other.atom))
    }
}

fn push_candidate(heap: &mut BinaryHeap<Candidate>, candidate: Candidate, count: usize) {
    heap.push(candidate);
    if heap.len() > count {
        let _discarded = heap.pop();
    }
}

fn sorted_candidates(heap: BinaryHeap<Candidate>) -> Vec<(u32, f32)> {
    let mut candidates: Vec<(u32, f32)> = heap
        .into_iter()
        .map(|candidate| (candidate.atom, candidate.squared))
        .collect();
    candidates.sort_unstable_by(|left, right| {
        left.1
            .total_cmp(&right.1)
            .then_with(|| left.0.cmp(&right.0))
    });
    candidates
}

fn nearest_brute(
    positions: &[[f32; 3]],
    targets: &[u32],
    query: u32,
    point: [f32; 3],
    count: usize,
    periodic: Option<&PeriodicBox>,
) -> Vec<(u32, f32)> {
    let mut candidates: Vec<(u32, f32)> = targets
        .iter()
        .filter(|target| **target != query)
        .filter_map(|target| {
            let position = positions.get(*target as usize).copied()?;
            finite(position).then(|| (*target, distance_squared(point, position, periodic)))
        })
        .collect();
    candidates.sort_unstable_by(|left, right| {
        left.1
            .total_cmp(&right.1)
            .then_with(|| left.0.cmp(&right.0))
    });
    candidates.truncate(count);
    candidates
}

#[cfg(test)]
#[path = "kd_tests.rs"]
mod tests;
