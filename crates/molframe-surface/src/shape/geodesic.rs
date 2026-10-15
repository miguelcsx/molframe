//! Deterministic shortest paths along surface-mesh edges.

use crate::{IndexedSurfaceMesh, numeric::usize_to_u32};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Why intrinsic traversal of a surface was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum SurfaceGeometryError {
    /// A source vertex is outside the mesh.
    #[error("surface vertex is out of bounds")]
    VertexOutOfBounds,
    /// The mesh has non-manifold adjacency or degenerate faces.
    #[error("surface mesh is non-manifold or degenerate")]
    NonManifoldMesh,
    /// A radius was negative or non-finite.
    #[error("surface radius is invalid")]
    InvalidRadius,
}

/// Distances from one source along mesh edges.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceDistances {
    /// Source vertex.
    pub source: u32,
    /// Per-vertex distance; infinity means a disconnected component.
    pub distances: Box<[f64]>,
}

/// Computes the graph-geodesic approximation from one source.
///
/// Runtime is `O((V + E) log V)` and output space is `O(V)`.
///
/// # Errors
///
/// Refuses invalid sources and non-manifold meshes.
pub fn edge_geodesic_distances(
    mesh: &IndexedSurfaceMesh,
    source: u32,
) -> Result<SurfaceDistances, SurfaceGeometryError> {
    if source as usize >= mesh.vertices.len() {
        return Err(SurfaceGeometryError::VertexOutOfBounds);
    }
    if !mesh.report.is_manifold() {
        return Err(SurfaceGeometryError::NonManifoldMesh);
    }
    let graph = weighted_graph(mesh);
    let mut distances = vec![f64::INFINITY; mesh.vertices.len()];
    distances[source as usize] = 0.0;
    let mut heap = BinaryHeap::from([QueueEntry {
        distance: 0.0,
        vertex: source,
    }]);
    while let Some(entry) = heap.pop() {
        if entry.distance > distances[entry.vertex as usize] {
            continue;
        }
        for &(next, weight) in &graph[entry.vertex as usize] {
            let candidate = entry.distance + weight;
            if candidate < distances[next as usize] {
                distances[next as usize] = candidate;
                heap.push(QueueEntry {
                    distance: candidate,
                    vertex: next,
                });
            }
        }
    }
    Ok(SurfaceDistances {
        source,
        distances: distances.into_boxed_slice(),
    })
}

/// Vertices within an intrinsic radius, ordered by vertex index.
///
/// # Errors
///
/// Refuses invalid radii or any error from geodesic traversal.
pub fn surface_patch(
    mesh: &IndexedSurfaceMesh,
    source: u32,
    radius: f64,
) -> Result<Box<[u32]>, SurfaceGeometryError> {
    if !radius.is_finite() || radius < 0.0 {
        return Err(SurfaceGeometryError::InvalidRadius);
    }
    let distances = edge_geodesic_distances(mesh, source)?;
    Ok(distances
        .distances
        .iter()
        .enumerate()
        .filter_map(|(index, distance)| {
            if *distance <= radius {
                Some(usize_to_u32(index))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
        .into_boxed_slice())
}

#[derive(Clone, Copy, Debug)]
struct QueueEntry {
    distance: f64,
    vertex: u32,
}

impl PartialEq for QueueEntry {
    fn eq(&self, other: &Self) -> bool {
        self.vertex == other.vertex && self.distance.to_bits() == other.distance.to_bits()
    }
}
impl Eq for QueueEntry {}
impl PartialOrd for QueueEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for QueueEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.vertex.cmp(&self.vertex))
    }
}

fn weighted_graph(mesh: &IndexedSurfaceMesh) -> Vec<Vec<(u32, f64)>> {
    let mut graph = vec![Vec::new(); mesh.vertices.len()];
    for face in &mesh.faces {
        for (left, right) in [
            (face.0[0], face.0[1]),
            (face.0[1], face.0[2]),
            (face.0[2], face.0[0]),
        ] {
            let Some(a) = mesh.vertices.get(left as usize) else {
                continue;
            };
            let Some(b) = mesh.vertices.get(right as usize) else {
                continue;
            };
            let weight = distance(*a, *b);
            graph[left as usize].push((right, weight));
            graph[right as usize].push((left, weight));
        }
    }
    for neighbours in &mut graph {
        neighbours.sort_by(|left, right| {
            left.0
                .cmp(&right.0)
                .then_with(|| left.1.total_cmp(&right.1))
        });
        neighbours.dedup_by(|left, right| left.0 == right.0);
    }
    graph
}

fn distance(left: [f32; 3], right: [f32; 3]) -> f64 {
    (0..3)
        .map(|axis| f64::from(right[axis] - left[axis]).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
#[path = "geodesic_tests.rs"]
mod tests;
