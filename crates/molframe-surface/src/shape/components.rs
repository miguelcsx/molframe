//! Deterministic connected components and presentation-oriented mesh filtering.
//!
//! Building vertex-to-face adjacency and walking every face costs `O(v + f)`
//! memory and time, excluding the final deterministic component sort.

use crate::{IndexedSurfaceMesh, SurfaceFace};
use std::collections::VecDeque;

#[cfg(test)]
#[path = "components_tests.rs"]
mod tests;

/// One face-connected surface component.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceComponent {
    /// Source face indices in ascending order.
    pub faces: Vec<u32>,
    /// Triangle area in squared coordinate units.
    pub area: f64,
}

/// Selection applied to connected surface components.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceComponentFilter {
    /// Components below this area are discarded.
    pub minimum_area: f64,
    /// Optional maximum number of largest components retained.
    pub maximum_components: Option<usize>,
}

impl Default for SurfaceComponentFilter {
    fn default() -> Self {
        Self {
            minimum_area: 0.0,
            maximum_components: None,
        }
    }
}

/// Why a surface-component filter could not be applied.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq)]
pub enum SurfaceComponentError {
    /// Area and count limits must be finite and meaningful.
    #[error("surface component limits must be finite and non-negative, with a positive count")]
    InvalidFilter,
    /// The compacted mesh exceeds the 32-bit index contract.
    #[error("filtered surface exceeds the 32-bit mesh index range")]
    MeshTooLarge,
}

/// Finds face-connected components in stable source-face order.
#[must_use]
pub fn surface_components(mesh: &IndexedSurfaceMesh) -> Vec<SurfaceComponent> {
    let mut incident = vec![Vec::<usize>::new(); mesh.vertices.len()];
    for (face_index, face) in mesh.faces.iter().enumerate() {
        for vertex in face.0 {
            if let Some(faces) = incident.get_mut(vertex as usize) {
                faces.push(face_index);
            }
        }
    }
    let mut visited = vec![false; mesh.faces.len()];
    let mut components = Vec::new();
    for start in 0..mesh.faces.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut queue = VecDeque::from([start]);
        let mut faces = Vec::new();
        let mut area = 0.0_f64;
        while let Some(face_index) = queue.pop_front() {
            let Some(face) = mesh.faces.get(face_index).copied() else {
                continue;
            };
            if let Ok(index) = u32::try_from(face_index) {
                faces.push(index);
            }
            area += face_area(mesh, face);
            for vertex in face.0 {
                let Some(neighbours) = incident.get(vertex as usize) else {
                    continue;
                };
                for neighbour in neighbours {
                    if !visited[*neighbour] {
                        visited[*neighbour] = true;
                        queue.push_back(*neighbour);
                    }
                }
            }
        }
        components.push(SurfaceComponent { faces, area });
    }
    components
}

/// Retains the largest components satisfying a minimum area and compacts them.
///
/// # Errors
///
/// Returns [`SurfaceComponentError`] for invalid limits or index overflow.
pub fn filter_surface_components(
    mesh: &IndexedSurfaceMesh,
    filter: SurfaceComponentFilter,
) -> Result<IndexedSurfaceMesh, SurfaceComponentError> {
    validate_filter(filter)?;
    let mut components = surface_components(mesh);
    components.retain(|component| component.area >= filter.minimum_area);
    components.sort_by(|left, right| {
        right
            .area
            .total_cmp(&left.area)
            .then_with(|| first_face(left).cmp(&first_face(right)))
    });
    if let Some(maximum) = filter.maximum_components {
        components.truncate(maximum);
    }
    let mut kept = components
        .iter()
        .flat_map(|component| component.faces.iter().copied())
        .collect::<Vec<_>>();
    kept.sort_unstable();
    compact_faces(mesh, &kept)
}

fn compact_faces(
    mesh: &IndexedSurfaceMesh,
    kept: &[u32],
) -> Result<IndexedSurfaceMesh, SurfaceComponentError> {
    let mut remap = vec![None; mesh.vertices.len()];
    let mut vertices = Vec::new();
    let mut faces = Vec::with_capacity(kept.len());
    for face_index in kept {
        let Some(source) = mesh.faces.get(*face_index as usize) else {
            continue;
        };
        let mut target = [0_u32; 3];
        let mut valid = true;
        for (corner, source_vertex) in source.0.into_iter().enumerate() {
            let Some(position) = mesh.vertices.get(source_vertex as usize).copied() else {
                valid = false;
                break;
            };
            let mapped = if let Some(index) = remap[source_vertex as usize] {
                index
            } else {
                let index = u32::try_from(vertices.len())
                    .map_err(|_| SurfaceComponentError::MeshTooLarge)?;
                vertices.push(position);
                remap[source_vertex as usize] = Some(index);
                index
            };
            target[corner] = mapped;
        }
        if valid {
            faces.push(SurfaceFace(target));
        }
    }
    Ok(IndexedSurfaceMesh::new(vertices, faces))
}

fn face_area(mesh: &IndexedSurfaceMesh, face: SurfaceFace) -> f64 {
    let Some(a) = mesh.vertices.get(face.0[0] as usize) else {
        return 0.0;
    };
    let Some(b) = mesh.vertices.get(face.0[1] as usize) else {
        return 0.0;
    };
    let Some(c) = mesh.vertices.get(face.0[2] as usize) else {
        return 0.0;
    };
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    cross
        .into_iter()
        .map(|value| f64::from(value) * f64::from(value))
        .sum::<f64>()
        .sqrt()
        * 0.5
}

fn validate_filter(filter: SurfaceComponentFilter) -> Result<(), SurfaceComponentError> {
    if !filter.minimum_area.is_finite()
        || filter.minimum_area < 0.0
        || matches!(filter.maximum_components, Some(0))
    {
        Err(SurfaceComponentError::InvalidFilter)
    } else {
        Ok(())
    }
}

fn first_face(component: &SurfaceComponent) -> u32 {
    let Some(face) = component.faces.first().copied() else {
        return u32::MAX;
    };
    face
}
