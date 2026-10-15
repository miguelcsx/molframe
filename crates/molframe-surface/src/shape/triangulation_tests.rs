use super::{IndexedSurfaceMesh, SurfaceFace};

#[test]
fn analytic_normals_preserve_the_provider_columns_instead_of_recomputing_face_normals() {
    let token = molframe_core::execution::CancellationToken::default();
    let vertices = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let faces = vec![SurfaceFace([0, 1, 2])];
    let normals = vec![[1.0, 0.0, 0.0]; 3];
    let pointer = normals.as_ptr();
    let geometry =
        IndexedSurfaceMesh::from_analytic_normals(vertices.clone(), faces.clone(), normals, &token)
            .expect("analytic geometry");
    assert_eq!(geometry.vertex_normals.as_ptr(), pointer);
    assert_eq!(
        geometry.report,
        IndexedSurfaceMesh::new(vertices.clone(), faces.clone()).report
    );
    assert!(
        IndexedSurfaceMesh::from_analytic_normals(
            vertices.clone(),
            faces.clone(),
            vec![[0.0; 3]; 3],
            &token,
        )
        .is_err()
    );
    assert!(
        IndexedSurfaceMesh::from_analytic_normals(
            vertices.clone(),
            faces.clone(),
            vec![[1.0, 0.0, 0.0]; 2],
            &token,
        )
        .is_err()
    );
    token.cancel();
    assert!(matches!(
        IndexedSurfaceMesh::from_analytic_normals(
            vertices,
            faces,
            vec![[1.0, 0.0, 0.0]; 3],
            &token,
        ),
        Err(crate::SurfaceOutputError::Cancelled(_))
    ));
}

#[test]
fn nonfinite_face_coordinates_are_reported_as_degenerate() {
    let geometry = IndexedSurfaceMesh::new(
        vec![[f32::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![SurfaceFace([0, 1, 2])],
    );
    assert_eq!(geometry.report.degenerate_faces, 1);
    assert!(!geometry.report.is_manifold());
}

#[test]
fn small_nonzero_triangles_keep_their_area_and_unit_normals() {
    let mesh = IndexedSurfaceMesh::new(
        vec![[0.0; 3], [1e-6, 0.0, 0.0], [0.0, 1e-6, 0.0]],
        vec![SurfaceFace([0, 1, 2])],
    );
    assert_eq!(mesh.report.degenerate_faces, 0);
    assert!((mesh.area() - 5e-13).abs() < 1e-19);
    for normal in &mesh.vertex_normals {
        assert_eq!(normal.map(f32::to_bits), [0.0, 0.0, 1.0].map(f32::to_bits));
    }
}

#[test]
fn vertex_connected_faces_share_one_component_and_unused_vertices_stay_separate() {
    let mesh = IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [10.0, 10.0, 10.0],
        ],
        vec![SurfaceFace([0, 1, 2]), SurfaceFace([0, 3, 4])],
    );
    assert_eq!(mesh.report.components, 2);
    assert_eq!(mesh.report.boundary_edges, 6);
    assert_eq!(mesh.report.degenerate_faces, 0);
}

#[test]
fn a_long_mesh_strip_has_one_component_without_recursive_traversal() {
    let vertices = (0..1000u16)
        .map(|row| [f32::from(row), f32::from(row % 2), 0.0])
        .collect();
    let faces = (0..998u32)
        .map(|row| SurfaceFace([row, row + 1, row + 2]))
        .collect();
    let mesh = IndexedSurfaceMesh::new(vertices, faces);
    assert_eq!(mesh.report.components, 1);
    assert_eq!(mesh.report.degenerate_faces, 0);
}

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

#[test]
fn finite_large_coordinates_keep_area_normals_and_degenerate_faces_well_defined() {
    let extent = f32::MAX;
    let triangle = IndexedSurfaceMesh::new(
        vec![[-extent, 0.0, 0.0], [extent, 0.0, 0.0], [0.0, extent, 0.0]],
        vec![SurfaceFace([0, 1, 2])],
    );
    assert_eq!(triangle.report.degenerate_faces, 0);
    assert!(triangle.area().is_finite() && triangle.area() > 0.0);
    assert!(
        triangle
            .vertex_normals
            .iter()
            .all(|normal| normal.map(f32::to_bits) == [0.0_f32, 0.0, 1.0].map(f32::to_bits))
    );
    let collinear = IndexedSurfaceMesh::new(
        vec![[-extent, 0.0, 0.0], [0.0, 0.0, 0.0], [extent, 0.0, 0.0]],
        vec![SurfaceFace([0, 1, 2])],
    );
    assert_eq!(collinear.report.degenerate_faces, 1);
    assert_eq!(collinear.area().to_bits(), 0.0_f64.to_bits());
    assert!(
        collinear
            .vertex_normals
            .iter()
            .all(|normal| normal.map(f32::to_bits) == [0.0_f32; 3].map(f32::to_bits))
    );
}
