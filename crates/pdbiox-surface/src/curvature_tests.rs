use super::{CurvatureQuality, surface_curvatures};
use crate::{IndexedSurfaceMesh, SurfaceFace};

#[test]
fn a_planar_interior_vertex_has_zero_curvature() {
    let mesh = IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        vec![
            SurfaceFace([0, 1, 2]),
            SurfaceFace([0, 2, 3]),
            SurfaceFace([0, 3, 4]),
            SurfaceFace([0, 4, 1]),
        ],
    );
    let values = surface_curvatures(&mesh);
    assert!(values.is_ok_and(|values| {
        values[0].quality == CurvatureQuality::Interior
            && values[0].mean.abs() < 1e-12
            && values[0].gaussian.abs() < 1e-12
    }));
}

#[test]
fn a_closed_octahedron_has_positive_gaussian_curvature() {
    let mesh = IndexedSurfaceMesh::new(
        vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
        vec![
            SurfaceFace([0, 2, 4]),
            SurfaceFace([2, 1, 4]),
            SurfaceFace([1, 3, 4]),
            SurfaceFace([3, 0, 4]),
            SurfaceFace([2, 0, 5]),
            SurfaceFace([1, 2, 5]),
            SurfaceFace([3, 1, 5]),
            SurfaceFace([0, 3, 5]),
        ],
    );
    let values = surface_curvatures(&mesh);
    assert!(values.is_ok_and(|values| values.iter().all(|value| value.gaussian > 0.0)));
}
