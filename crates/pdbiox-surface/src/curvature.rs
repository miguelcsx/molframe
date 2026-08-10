//! Discrete curvature on an indexed triangular surface.
//!
//! Gaussian curvature uses angle deficit and mean curvature uses the cotangent
//! Laplacian. Both are intrinsic to the mesh approximation, so callers must
//! retain mesh resolution and per-vertex quality when interpreting the values.

use crate::{IndexedSurfaceMesh, SurfaceFace, SurfaceGeometryError, triangulation::edge_incidence};
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, PI};

/// Reliability classification for one curvature estimate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurvatureQuality {
    /// Closed, manifold neighborhood with finite area.
    Interior,
    /// The neighborhood meets a mesh boundary.
    Boundary,
    /// The neighborhood is degenerate or cannot support an estimate.
    Degenerate,
}

/// Curvature descriptors at one surface vertex.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceCurvature {
    /// Mean curvature magnitude in Å⁻¹.
    pub mean: f64,
    /// Gaussian curvature in Å⁻².
    pub gaussian: f64,
    /// Larger principal curvature in Å⁻¹.
    pub maximum: f64,
    /// Smaller principal curvature in Å⁻¹.
    pub minimum: f64,
    /// Koenderink shape index in `[-1, 1]`.
    pub shape_index: f64,
    /// Root-mean-square principal curvature in Å⁻¹.
    pub curvedness: f64,
    /// Neighborhood quality.
    pub quality: CurvatureQuality,
}

/// Computes one curvature descriptor per mesh vertex.
///
/// # Errors
///
/// Refuses non-manifold or degenerate mesh topology.
pub fn surface_curvatures(
    mesh: &IndexedSurfaceMesh,
) -> Result<Box<[SurfaceCurvature]>, SurfaceGeometryError> {
    if !mesh.report.is_manifold() {
        return Err(SurfaceGeometryError::NonManifoldMesh);
    }
    let edges = edge_incidence(&mesh.faces);
    let mut areas = vec![0.0; mesh.vertices.len()];
    let mut angle_sums = vec![0.0; mesh.vertices.len()];
    let mut cotangent_weights = BTreeMap::new();
    for face in &mesh.faces {
        accumulate_face(
            mesh,
            *face,
            &mut areas,
            &mut angle_sums,
            &mut cotangent_weights,
        );
    }
    Ok((0..mesh.vertices.len())
        .map(|vertex| {
            curvature_at(
                mesh,
                vertex,
                areas[vertex],
                angle_sums[vertex],
                &edges,
                &cotangent_weights,
            )
        })
        .collect::<Vec<_>>()
        .into_boxed_slice())
}

fn accumulate_face(
    mesh: &IndexedSurfaceMesh,
    face: SurfaceFace,
    areas: &mut [f64],
    angle_sums: &mut [f64],
    weights: &mut BTreeMap<(u32, u32), f64>,
) {
    let (Some(a), Some(b), Some(c)) = (
        mesh.vertices.get(face.0[0] as usize),
        mesh.vertices.get(face.0[1] as usize),
        mesh.vertices.get(face.0[2] as usize),
    ) else {
        return;
    };
    let points = [to_f64(*a), to_f64(*b), to_f64(*c)];
    let area = 0.5
        * length(cross(
            subtract(points[1], points[0]),
            subtract(points[2], points[0]),
        ));
    if !area.is_finite() || area <= f64::EPSILON {
        return;
    }
    for local in 0..3 {
        let vertex = face.0[local] as usize;
        areas[vertex] += area / 3.0;
        let left = subtract(points[(local + 1) % 3], points[local]);
        let right = subtract(points[(local + 2) % 3], points[local]);
        let angle = angle_between(left, right);
        angle_sums[vertex] += angle;
        let opposite = (face.0[(local + 1) % 3], face.0[(local + 2) % 3]);
        let edge = (opposite.0.min(opposite.1), opposite.0.max(opposite.1));
        let cotangent = angle.cos() / angle.sin().max(f64::EPSILON);
        let current = match weights.get(&edge).copied() {
            Some(value) => value,
            None => 0.0,
        };
        weights.insert(edge, current + cotangent);
    }
}

fn curvature_at(
    mesh: &IndexedSurfaceMesh,
    vertex: usize,
    area: f64,
    angle_sum: f64,
    edges: &BTreeMap<(u32, u32), u32>,
    weights: &BTreeMap<(u32, u32), f64>,
) -> SurfaceCurvature {
    let boundary = edges.iter().any(|(&(left, right), &count)| {
        count == 1 && (left as usize == vertex || right as usize == vertex)
    });
    if area <= f64::EPSILON || !area.is_finite() {
        return zero(CurvatureQuality::Degenerate);
    }
    let mut laplacian = [0.0; 3];
    for (&(left, right), &weight) in weights {
        let neighbour = if left as usize == vertex {
            Some(right)
        } else if right as usize == vertex {
            Some(left)
        } else {
            None
        };
        let Some(neighbour) = neighbour else {
            continue;
        };
        let Some(origin) = mesh.vertices.get(vertex) else {
            continue;
        };
        let Some(target) = mesh.vertices.get(neighbour as usize) else {
            continue;
        };
        for axis in 0..3 {
            laplacian[axis] += weight * f64::from(target[axis] - origin[axis]);
        }
    }
    let mean = length(laplacian) / (4.0 * area);
    let target_angle = if boundary { PI } else { 2.0 * PI };
    let gaussian = (target_angle - angle_sum) / area;
    let discriminant = mean.mul_add(mean, -gaussian).max(0.0).sqrt();
    let maximum = mean + discriminant;
    let minimum = mean - discriminant;
    let difference = maximum - minimum;
    let ratio = if difference.abs() <= f64::EPSILON {
        if maximum + minimum >= 0.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        }
    } else {
        (maximum + minimum) / difference
    };
    SurfaceCurvature {
        mean,
        gaussian,
        maximum,
        minimum,
        shape_index: 2.0 / PI * ratio.atan(),
        curvedness: ((maximum * maximum + minimum * minimum) * 0.5).sqrt(),
        quality: if boundary {
            CurvatureQuality::Boundary
        } else {
            CurvatureQuality::Interior
        },
    }
}

fn zero(quality: CurvatureQuality) -> SurfaceCurvature {
    SurfaceCurvature {
        mean: 0.0,
        gaussian: 0.0,
        maximum: 0.0,
        minimum: 0.0,
        shape_index: 0.0,
        curvedness: 0.0,
        quality,
    }
}

fn to_f64(point: [f32; 3]) -> [f64; 3] {
    point.map(f64::from)
}
fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}
fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}
fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}
fn length(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}
fn angle_between(left: [f64; 3], right: [f64; 3]) -> f64 {
    let denominator = length(left) * length(right);
    if denominator <= f64::EPSILON {
        return FRAC_PI_2;
    }
    (dot(left, right) / denominator).clamp(-1.0, 1.0).acos()
}

#[cfg(test)]
#[path = "curvature_tests.rs"]
mod tests;
