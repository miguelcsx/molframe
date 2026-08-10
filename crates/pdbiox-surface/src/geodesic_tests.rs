use super::{edge_geodesic_distances, surface_patch};
use crate::{IndexedSurfaceMesh, SurfaceFace};

fn square() -> IndexedSurfaceMesh {
    IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![SurfaceFace([0, 1, 2]), SurfaceFace([0, 2, 3])],
    )
}

#[test]
fn edge_distance_follows_the_shortest_mesh_path() {
    let distances = edge_geodesic_distances(&square(), 1);
    assert!(distances.is_ok_and(|value| (value.distances[3] - 2.0).abs() < 1e-12));
}

#[test]
fn a_patch_is_sorted_and_radius_bounded() {
    let patch = surface_patch(&square(), 0, 1.01);
    assert_eq!(patch.as_deref(), Ok(&[0, 1, 3][..]));
}
