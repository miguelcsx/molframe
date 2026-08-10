//! Deterministic clustering and representatives over ensemble distance matrices.

use crate::numeric::f64_from_usize;
use crate::{EnsembleDistanceMatrix, EnsembleGeometryError};
use std::collections::VecDeque;

/// Agglomerative inter-cluster distance definition.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Linkage {
    /// Minimum pair distance.
    Single,
    /// Maximum pair distance.
    Complete,
    /// Arithmetic mean pair distance.
    Average,
}

/// Deterministic cluster membership and medoid representatives.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Clustering {
    /// Cluster label per observation; `None` denotes DBSCAN noise.
    pub labels: Vec<Option<usize>>,
    /// Sorted observation indices for each cluster.
    pub members: Vec<Vec<usize>>,
    /// Minimum-total-distance member of each cluster.
    pub medoids: Vec<usize>,
}

/// Agglomerative clustering to an explicit number of clusters.
///
/// Ties are resolved by the lexicographically smallest pair of current member
/// lists, making results independent of hash iteration or thread scheduling.
///
/// # Errors
///
/// Returns [`EnsembleGeometryError`] for an invalid matrix or cluster count.
pub fn agglomerative_clustering(
    distances: &EnsembleDistanceMatrix,
    cluster_count: usize,
    linkage: Linkage,
) -> Result<Clustering, EnsembleGeometryError> {
    validate_matrix(distances)?;
    if cluster_count == 0 || cluster_count > distances.size {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let mut clusters: Vec<Vec<usize>> = (0..distances.size).map(|index| vec![index]).collect();
    while clusters.len() > cluster_count {
        let mut best: Option<(f64, usize, usize)> = None;
        for left in 0..clusters.len() {
            for right in left + 1..clusters.len() {
                let distance =
                    linkage_distance(distances, &clusters[left], &clusters[right], linkage)?;
                if best.is_none_or(|candidate| (distance, left, right) < candidate) {
                    best = Some((distance, left, right));
                }
            }
        }
        let Some((_, left, right)) = best else {
            return Err(EnsembleGeometryError::InvalidParameter);
        };
        let removed = clusters.remove(right);
        clusters[left].extend(removed);
        clusters[left].sort_unstable();
    }
    clusters.sort();
    make_clustering(distances, clusters, &[])
}

/// Density-based clustering with explicit neighbourhood and core-size policy.
///
/// # Errors
///
/// Returns [`EnsembleGeometryError`] for an invalid matrix, epsilon, or core size.
pub fn dbscan_clustering(
    distances: &EnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
) -> Result<Clustering, EnsembleGeometryError> {
    validate_matrix(distances)?;
    if !epsilon.is_finite() || epsilon < 0.0 || minimum_points == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let mut labels = vec![None; distances.size];
    let mut visited = vec![false; distances.size];
    let mut cluster_count = 0;
    for point in 0..distances.size {
        if visited[point] {
            continue;
        }
        visited[point] = true;
        let neighbours = neighbourhood(distances, point, epsilon);
        if neighbours.len() < minimum_points {
            continue;
        }
        let mut expansion = ClusterExpansion {
            distances,
            epsilon,
            minimum_points,
            label: cluster_count,
            visited: &mut visited,
            labels: &mut labels,
        };
        expansion.run(neighbours);
        cluster_count += 1;
    }
    let mut members = vec![Vec::new(); cluster_count];
    for (point, label) in labels.iter().copied().enumerate() {
        if let Some(label) = label {
            members[label].push(point);
        }
    }
    let medoids = members
        .iter()
        .map(|cluster| medoid(distances, cluster))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Clustering {
        labels,
        members,
        medoids,
    })
}

/// Returns the member minimizing total within-group distance.
///
/// # Errors
///
/// Returns [`EnsembleGeometryError`] for invalid matrix or empty/out-of-range members.
pub fn medoid(
    distances: &EnsembleDistanceMatrix,
    members: &[usize],
) -> Result<usize, EnsembleGeometryError> {
    validate_matrix(distances)?;
    if members.is_empty() || members.iter().any(|member| *member >= distances.size) {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    members
        .iter()
        .copied()
        .map(|candidate| {
            let total = members
                .iter()
                .map(|other| distances.values[candidate * distances.size + other])
                .sum::<f64>();
            (total, candidate)
        })
        .min_by(|left, right| left.0.total_cmp(&right.0).then(left.1.cmp(&right.1)))
        .map(|(_, candidate)| candidate)
        .ok_or(EnsembleGeometryError::InvalidParameter)
}

fn validate_matrix(distances: &EnsembleDistanceMatrix) -> Result<(), EnsembleGeometryError> {
    if distances.size == 0 || distances.values.len() != distances.size * distances.size {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }
    for row in 0..distances.size {
        for column in 0..distances.size {
            let value = distances.values[row * distances.size + column];
            let reverse = distances.values[column * distances.size + row];
            if !value.is_finite() || value < 0.0 || value.to_bits() != reverse.to_bits() {
                return Err(EnsembleGeometryError::InvalidParameter);
            }
        }
    }
    Ok(())
}

fn linkage_distance(
    distances: &EnsembleDistanceMatrix,
    left: &[usize],
    right: &[usize],
    linkage: Linkage,
) -> Result<f64, EnsembleGeometryError> {
    let values = left.iter().flat_map(|first| {
        right
            .iter()
            .map(move |second| distances.values[first * distances.size + second])
    });
    let distance = match linkage {
        Linkage::Single => values.fold(f64::INFINITY, f64::min),
        Linkage::Complete => values.fold(0.0, f64::max),
        Linkage::Average => {
            let pairs = left
                .len()
                .checked_mul(right.len())
                .and_then(f64_from_usize)
                .ok_or(EnsembleGeometryError::InvalidParameter)?;
            values.sum::<f64>() / pairs
        }
    };
    Ok(distance)
}

fn make_clustering(
    distances: &EnsembleDistanceMatrix,
    members: Vec<Vec<usize>>,
    noise: &[usize],
) -> Result<Clustering, EnsembleGeometryError> {
    let mut labels = vec![None; distances.size];
    for (label, cluster) in members.iter().enumerate() {
        for &point in cluster {
            labels[point] = Some(label);
        }
    }
    for &point in noise {
        labels[point] = None;
    }
    let medoids = members
        .iter()
        .map(|cluster| medoid(distances, cluster))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Clustering {
        labels,
        members,
        medoids,
    })
}

fn neighbourhood(distances: &EnsembleDistanceMatrix, point: usize, epsilon: f64) -> Vec<usize> {
    (0..distances.size)
        .filter(|other| distances.values[point * distances.size + other] <= epsilon)
        .collect()
}

struct ClusterExpansion<'a> {
    distances: &'a EnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
    label: usize,
    visited: &'a mut [bool],
    labels: &'a mut [Option<usize>],
}

impl ClusterExpansion<'_> {
    fn run(&mut self, neighbours: Vec<usize>) {
        let mut queue: VecDeque<usize> = neighbours.into();
        while let Some(point) = queue.pop_front() {
            if !self.visited[point] {
                self.visited[point] = true;
                let neighbours = neighbourhood(self.distances, point, self.epsilon);
                if neighbours.len() >= self.minimum_points {
                    for neighbour in neighbours {
                        if !queue.contains(&neighbour) {
                            queue.push_back(neighbour);
                        }
                    }
                }
            }
            self.labels[point].get_or_insert(self.label);
        }
    }
}

#[cfg(test)]
#[path = "clustering_tests.rs"]
mod tests;
