//! Density-based clustering over a distance matrix.
//!
//! Complexity: O(n^2) time from neighbourhood scans over n frames; O(n) extra space.

use super::matrix::{matrix_row, reusable_index_buffer, validate_matrix};
use super::medoid::cluster_medoids;
use super::types::Clustering;
use crate::{EnsembleDistanceMatrix, EnsembleGeometryError};

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

    let size = distances.size;
    let mut labels = vec![None; size];
    let mut visited = vec![false; size];
    let mut queued_epoch = vec![usize::MAX; size];

    let mut neighbours = reusable_index_buffer(size)?;
    let mut queue = reusable_index_buffer(size)?;
    let mut scratch = reusable_index_buffer(size)?;

    let mut cluster_count = 0usize;

    for point in 0..size {
        let visited_slot = visited
            .get_mut(point)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        if *visited_slot {
            continue;
        }

        *visited_slot = true;
        collect_neighbourhood(distances, point, epsilon, &mut neighbours)?;

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
            queued_epoch: &mut queued_epoch,
        };

        expansion.run(&neighbours, &mut queue, &mut scratch)?;

        cluster_count = cluster_count
            .checked_add(1)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
    }

    let mut members = vec![Vec::new(); cluster_count];

    for (point, label) in labels.iter().copied().enumerate() {
        let Some(label) = label else {
            continue;
        };

        let cluster = members
            .get_mut(label)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;
        cluster.push(point);
    }

    let medoids = cluster_medoids(distances, &members)?;

    Ok(Clustering {
        labels,
        members,
        medoids,
    })
}

/// Collects every observation within `epsilon` of `point` into reusable storage.
fn collect_neighbourhood(
    distances: &EnsembleDistanceMatrix,
    point: usize,
    epsilon: f64,
    output: &mut Vec<usize>,
) -> Result<(), EnsembleGeometryError> {
    output.clear();

    let row = matrix_row(distances, point)?;

    for (other, &distance) in row.iter().enumerate() {
        if distance <= epsilon {
            output.push(other);
        }
    }

    Ok(())
}

/// State required to expand one DBSCAN cluster.
struct ClusterExpansion<'a> {
    distances: &'a EnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
    label: usize,
    visited: &'a mut [bool],
    labels: &'a mut [Option<usize>],
    queued_epoch: &'a mut [usize],
}

impl ClusterExpansion<'_> {
    /// Expands one core neighbourhood and assigns previously unlabeled points.
    fn run(
        &mut self,
        neighbours: &[usize],
        queue: &mut Vec<usize>,
        scratch: &mut Vec<usize>,
    ) -> Result<(), EnsembleGeometryError> {
        queue.clear();

        for &point in neighbours {
            self.enqueue(point, queue)?;
        }

        let mut head = 0usize;

        while let Some(&point) = queue.get(head) {
            head += 1;

            let was_visited = self
                .visited
                .get(point)
                .copied()
                .ok_or(EnsembleGeometryError::InvalidParameter)?;

            if !was_visited {
                let visited = self
                    .visited
                    .get_mut(point)
                    .ok_or(EnsembleGeometryError::InvalidParameter)?;
                *visited = true;

                collect_neighbourhood(self.distances, point, self.epsilon, scratch)?;

                if scratch.len() >= self.minimum_points {
                    for &neighbour in scratch.iter() {
                        self.enqueue(neighbour, queue)?;
                    }
                }
            }

            let label = self
                .labels
                .get_mut(point)
                .ok_or(EnsembleGeometryError::InvalidParameter)?;

            if label.is_none() {
                *label = Some(self.label);
            }
        }

        Ok(())
    }

    /// Adds a point once to the current cluster's expansion queue.
    fn enqueue(
        &mut self,
        point: usize,
        queue: &mut Vec<usize>,
    ) -> Result<(), EnsembleGeometryError> {
        let epoch = self
            .queued_epoch
            .get_mut(point)
            .ok_or(EnsembleGeometryError::InvalidParameter)?;

        if *epoch == self.label {
            return Ok(());
        }

        *epoch = self.label;
        queue.push(point);
        Ok(())
    }
}
