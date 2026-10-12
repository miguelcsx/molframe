//! Flat deterministic edge incidence and disjoint-set topology diagnostics.
use crate::{
    MeshReport, SurfaceFace,
    numeric::usize_to_u32,
    triangulation::{component_root, face_degenerate},
};

pub(crate) struct MeshDiagnostics {
    edges: Vec<(u32, u32)>,
    parents: Vec<usize>,
    ranks: Vec<u8>,
    components: usize,
    degenerate_faces: usize,
}

impl MeshDiagnostics {
    pub(crate) fn new(vertices: usize, faces: usize) -> Self {
        Self {
            edges: Vec::with_capacity(faces.saturating_mul(3)),
            parents: (0..vertices).collect(),
            ranks: vec![0; vertices],
            components: vertices,
            degenerate_faces: 0,
        }
    }

    pub(crate) fn add(&mut self, vertices: &[[f32; 3]], face: SurfaceFace) {
        self.degenerate_faces += usize::from(face_degenerate(vertices, face));
        for (left, right) in [
            (face.0[0], face.0[1]),
            (face.0[1], face.0[2]),
            (face.0[2], face.0[0]),
        ] {
            self.edges.push((left.min(right), left.max(right)));
            self.join(left as usize, right as usize);
        }
    }

    fn join(&mut self, left: usize, right: usize) {
        if left >= self.parents.len() || right >= self.parents.len() {
            return;
        }
        let mut left = component_root(&mut self.parents, left);
        let mut right = component_root(&mut self.parents, right);
        if left == right {
            return;
        }
        if self.ranks[left] < self.ranks[right]
            || (self.ranks[left] == self.ranks[right] && left > right)
        {
            std::mem::swap(&mut left, &mut right);
        }
        self.parents[right] = left;
        if self.ranks[left] == self.ranks[right] {
            self.ranks[left] += 1;
        }
        self.components -= 1;
    }

    pub(crate) fn finish(mut self) -> MeshReport {
        // Sorting a contiguous edge column costs O(faces log faces), without
        // a tree-node allocation per edge or an adjacency vector per vertex.
        self.edges.sort_unstable();
        let mut boundary = 0;
        let mut non_manifold = 0;
        let mut first = 0;
        while first < self.edges.len() {
            let mut end = first + 1;
            while end < self.edges.len() && self.edges[end] == self.edges[first] {
                end += 1;
            }
            boundary += usize::from(end - first == 1);
            non_manifold += usize::from(end - first > 2);
            first = end;
        }
        MeshReport {
            boundary_edges: usize_to_u32(boundary),
            non_manifold_edges: usize_to_u32(non_manifold),
            degenerate_faces: usize_to_u32(self.degenerate_faces),
            components: usize_to_u32(self.components),
        }
    }
}
