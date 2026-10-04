use super::{IndexedSurfaceMesh, SurfaceFace};

#[test]
fn two_triangles_weld_their_shared_vertices() {
    let mesh = IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![SurfaceFace([0, 1, 2]), SurfaceFace([0, 2, 3])],
    );
    assert_eq!(mesh.report.boundary_edges, 4);
    assert_eq!(mesh.report.non_manifold_edges, 0);
    assert_eq!(mesh.report.components, 1);
    assert!(mesh.vertex_normals.iter().all(|normal| normal[2] > 0.99));
}

#[test]
fn an_edge_with_three_faces_is_non_manifold() {
    let mesh = IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        vec![
            SurfaceFace([0, 1, 2]),
            SurfaceFace([1, 0, 3]),
            SurfaceFace([0, 1, 4]),
        ],
    );
    assert_eq!(mesh.report.non_manifold_edges, 1);
    assert!(!mesh.report.is_manifold());
}

#[test]
fn the_area_of_a_unit_square_of_two_triangles_is_one() {
    let mesh = IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![SurfaceFace([0, 1, 2]), SurfaceFace([0, 2, 3])],
    );
    assert!((mesh.area() - 1.0).abs() < 1e-12);
}

#[test]
fn a_face_naming_a_missing_vertex_has_no_area_and_spoils_manifoldness() {
    let mesh = IndexedSurfaceMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![SurfaceFace([0, 1, 2]), SurfaceFace([0, 1, 9])],
    );
    assert!((mesh.area() - 0.5).abs() < 1e-12);
    assert!(!mesh.report.is_manifold());
}
