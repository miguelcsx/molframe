//! First-class parallel nearest-neighbour queries.

use super::indices::IndexWorkspace;
use crate::brute::distance_squared;
use crate::{
    KdPeriodicOptions, KdTree, KnnOptions, KnnResult, NO_NEIGHBOR, PeriodicBox, SpatialBackend,
    SpatialError,
};
use molframe_core::{AtomSelection, ExecutionContext};

/// Finds the `k` nearest target atoms of every query atom.
///
/// Rows follow the query selection's order. Within a row neighbours increase in
/// distance, equal distances ordering by atom index, and unfilled slots hold
/// [`NO_NEIGHBOR`] with an infinite squared distance. Periodic queries report
/// exact minimum-image distances; an unbounded periodic search looks through
/// two lattice shells of images, which covers every wrapped target of a cell
/// that is not extremely skewed, and a `max_distance` widens it exactly.
///
/// `BruteForce` scans every target serially and exists as a reference; every
/// other backend uses the parallel k-d tree, since a cell list has no
/// radius-free nearest query. Distances are f32; a target whose f64 distance
/// ties the k-th at the f32 rounding boundary may be ranked either way.
///
/// # Errors
///
/// Returns an error for an invalid `max_distance`, out-of-range atom indices,
/// memory admission failure or cancellation.
pub fn knn_within(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    targets: &AtomSelection,
    k: usize,
    options: KnnOptions,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<KnnResult, SpatialError> {
    let indices = IndexWorkspace::new(query, targets, positions.len(), context)?;
    if backend == SpatialBackend::BruteForce {
        return brute_knn(
            positions,
            indices.left(),
            indices.right(),
            k,
            options,
            periodic,
        );
    }
    KdTree::build_in(
        positions,
        indices.right(),
        periodic.copied(),
        KdPeriodicOptions::default(),
        context,
    )?
    .knn(indices.left(), k, options, context)
}

fn brute_knn(
    positions: &[[f32; 3]],
    query: &[u32],
    targets: &[u32],
    k: usize,
    options: KnnOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<KnnResult, SpatialError> {
    if let Some(distance) = options.max_distance
        && !(distance.is_finite() && distance >= 0.0)
    {
        return Err(SpatialError::InvalidCutoff);
    }
    let limit = options.max_distance.map_or(f32::INFINITY, |d| d * d);
    let mut indices = Vec::with_capacity(query.len() * k);
    let mut distances = Vec::with_capacity(query.len() * k);
    let finite = |point: [f32; 3]| point.iter().all(|value| value.is_finite());
    let mut row: Vec<(f32, u32)> = Vec::new();
    for &atom in query {
        let point = positions
            .get(atom as usize)
            .copied()
            .ok_or(SpatialError::AtomOutOfBounds(atom))?;
        row.clear();
        if finite(point) {
            for &target in targets {
                let other = positions
                    .get(target as usize)
                    .copied()
                    .ok_or(SpatialError::AtomOutOfBounds(target))?;
                if (options.exclude_self && target == atom) || !finite(other) {
                    continue;
                }
                let squared = distance_squared(point, other, periodic);
                if squared <= limit {
                    row.push((squared, target));
                }
            }
        }
        row.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        for slot in 0..k {
            let (id, squared) = row
                .get(slot)
                .map_or((NO_NEIGHBOR, f32::INFINITY), |&(d, id)| (id, d));
            indices.push(id);
            distances.push(squared);
        }
    }
    Ok(KnnResult {
        k,
        indices,
        distances_squared: distances,
    })
}
