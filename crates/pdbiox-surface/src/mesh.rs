//! Indexed topology derived from deterministic surface triangles.
//!
//! The SES writer keeps triangles independently because that is the natural
//! polygonisation output. Intrinsic algorithms instead need shared vertices and
//! edges. Conversion welds bit-identical vertices in first-seen order, validates
//! adjacency, and never changes the original surface representation.

use crate::{
    SolventExcludedSurface, SurfaceTriangle,
    numeric::{f64_to_f32, usize_to_u32},
};
use std::collections::{BTreeMap, VecDeque};

/// One indexed triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SurfaceFace(pub [u32; 3]);

/// Topological quality of an indexed surface.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MeshReport {
    /// Edges incident to exactly one face.
    pub boundary_edges: u32,
    /// Edges incident to more than two faces.
    pub non_manifold_edges: u32,
    /// Faces with repeated vertices or negligible area.
    pub degenerate_faces: u32,
    /// Connected face/vertex components.
    pub components: u32,
}

impl MeshReport {
    /// Whether every edge has at most two incident faces and all faces are valid.
    #[must_use]
    pub const fn is_manifold(&self) -> bool {
        self.non_manifold_edges == 0 && self.degenerate_faces == 0
    }
}

/// Shared-vertex triangular surface mesh.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IndexedSurfaceMesh {
    /// Vertex coordinates in first-seen triangle order.
    pub vertices: Vec<[f32; 3]>,
    /// Indexed triangular faces.
    pub faces: Vec<SurfaceFace>,
    /// Area-weighted unit vertex normals.
    pub vertex_normals: Vec<[f32; 3]>,
    /// Mesh topology diagnostics.
    pub report: MeshReport,
}

impl IndexedSurfaceMesh {
    /// Builds a shared-vertex mesh from independent triangles.
    #[must_use]
    pub fn from_triangles(triangles: &[SurfaceTriangle]) -> Self {
        let mut vertices = Vec::new();
        let mut lookup = BTreeMap::new();
        let mut faces = Vec::with_capacity(triangles.len());
        for triangle in triangles {
            let indices = triangle.vertices.map(|vertex| {
                let key = vertex.map(f32::to_bits);
                if let Some(index) = lookup.get(&key).copied() {
                    index
                } else {
                    let index = usize_to_u32(vertices.len());
                    vertices.push(vertex);
                    lookup.insert(key, index);
                    index
                }
            });
            faces.push(SurfaceFace(indices));
        }
        Self::new(vertices, faces)
    }

    /// Validates caller-supplied indexed geometry and computes normals.
    #[must_use]
    pub fn new(vertices: Vec<[f32; 3]>, faces: Vec<SurfaceFace>) -> Self {
        let adjacency = edge_incidence(&faces);
        let report = mesh_report(&vertices, &faces, &adjacency);
        let vertex_normals = vertex_normals(&vertices, &faces);
        Self {
            vertices,
            faces,
            vertex_normals,
            report,
        }
    }

    /// Sorted unique neighbours of one vertex.
    #[must_use]
    pub fn neighbours(&self, vertex: u32) -> Vec<u32> {
        let mut neighbours = Vec::new();
        for face in &self.faces {
            if face.0.contains(&vertex) {
                for candidate in face.0 {
                    if candidate != vertex {
                        neighbours.push(candidate);
                    }
                }
            }
        }
        neighbours.sort_unstable();
        neighbours.dedup();
        neighbours
    }

    /// Whether a vertex touches a boundary edge.
    #[must_use]
    pub fn is_boundary_vertex(&self, vertex: u32) -> bool {
        edge_incidence(&self.faces)
            .into_iter()
            .any(|((left, right), count)| count == 1 && (left == vertex || right == vertex))
    }
}

impl SolventExcludedSurface {
    /// Converts this surface to an indexed topology without changing it.
    #[must_use]
    pub fn indexed_mesh(&self) -> IndexedSurfaceMesh {
        IndexedSurfaceMesh::from_triangles(&self.triangles)
    }
}

pub(crate) fn edge_incidence(faces: &[SurfaceFace]) -> BTreeMap<(u32, u32), u32> {
    let mut edges = BTreeMap::new();
    for face in faces {
        for (left, right) in [
            (face.0[0], face.0[1]),
            (face.0[1], face.0[2]),
            (face.0[2], face.0[0]),
        ] {
            let edge = (left.min(right), left.max(right));
            let count = match edges.get(&edge).copied() {
                Some(value) => value,
                None => 0,
            };
            edges.insert(edge, count + 1);
        }
    }
    edges
}

fn mesh_report(
    vertices: &[[f32; 3]],
    faces: &[SurfaceFace],
    edges: &BTreeMap<(u32, u32), u32>,
) -> MeshReport {
    let boundary_edges = usize_to_u32(edges.values().filter(|count| **count == 1).count());
    let non_manifold_edges = usize_to_u32(edges.values().filter(|count| **count > 2).count());
    let degenerate_faces = usize_to_u32(
        faces
            .iter()
            .filter(|face| face_degenerate(vertices, **face))
            .count(),
    );
    let components = usize_to_u32(component_count(vertices.len(), edges));
    MeshReport {
        boundary_edges,
        non_manifold_edges,
        degenerate_faces,
        components,
    }
}

fn face_degenerate(vertices: &[[f32; 3]], face: SurfaceFace) -> bool {
    if face.0[0] == face.0[1] || face.0[1] == face.0[2] || face.0[2] == face.0[0] {
        return true;
    }
    let Some(a) = vertices.get(face.0[0] as usize) else {
        return true;
    };
    let Some(b) = vertices.get(face.0[1] as usize) else {
        return true;
    };
    let Some(c) = vertices.get(face.0[2] as usize) else {
        return true;
    };
    triangle_cross(*a, *b, *c)
        .iter()
        .map(|value| value * value)
        .sum::<f64>()
        <= f64::EPSILON
}

fn component_count(vertex_count: usize, edges: &BTreeMap<(u32, u32), u32>) -> usize {
    if vertex_count == 0 {
        return 0;
    }
    let mut graph = vec![Vec::new(); vertex_count];
    for &(left, right) in edges.keys() {
        if let (Some(left_edges), Some(_)) = (graph.get(left as usize), graph.get(right as usize)) {
            let _ = left_edges;
            graph[left as usize].push(right);
            graph[right as usize].push(left);
        }
    }
    let mut visited = vec![false; vertex_count];
    let mut components = 0;
    for start in 0..vertex_count {
        if visited[start] {
            continue;
        }
        components += 1;
        visited[start] = true;
        let mut queue = VecDeque::from([usize_to_u32(start)]);
        while let Some(vertex) = queue.pop_front() {
            for &next in &graph[vertex as usize] {
                if !visited[next as usize] {
                    visited[next as usize] = true;
                    queue.push_back(next);
                }
            }
        }
    }
    components
}

fn vertex_normals(vertices: &[[f32; 3]], faces: &[SurfaceFace]) -> Vec<[f32; 3]> {
    let mut normals = vec![[0.0_f64; 3]; vertices.len()];
    for face in faces {
        let (Some(a), Some(b), Some(c)) = (
            vertices.get(face.0[0] as usize),
            vertices.get(face.0[1] as usize),
            vertices.get(face.0[2] as usize),
        ) else {
            continue;
        };
        let cross = triangle_cross(*a, *b, *c);
        for vertex in face.0 {
            if let Some(normal) = normals.get_mut(vertex as usize) {
                for axis in 0..3 {
                    normal[axis] += cross[axis];
                }
            }
        }
    }
    normals
        .into_iter()
        .map(|normal| {
            let length = normal.iter().map(|value| value * value).sum::<f64>().sqrt();
            if length <= f64::EPSILON {
                [0.0; 3]
            } else {
                normal.map(|value| f64_to_f32(value / length))
            }
        })
        .collect()
}

pub(crate) fn triangle_cross(a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> [f64; 3] {
    let ab = [
        f64::from(b[0] - a[0]),
        f64::from(b[1] - a[1]),
        f64::from(b[2] - a[2]),
    ];
    let ac = [
        f64::from(c[0] - a[0]),
        f64::from(c[1] - a[1]),
        f64::from(c[2] - a[2]),
    ];
    [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ]
}

#[cfg(test)]
#[path = "mesh_tests.rs"]
mod tests;
