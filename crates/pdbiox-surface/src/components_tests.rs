use super::{
    SurfaceComponentError, SurfaceComponentFilter, filter_surface_components, surface_components,
};
use crate::{IndexedSurfaceMesh, SurfaceFace};

fn two_triangles() -> IndexedSurfaceMesh {
    IndexedSurfaceMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.5, 0.0, 0.0],
            [10.0, 0.5, 0.0],
        ],
        vec![SurfaceFace([0, 1, 2]), SurfaceFace([3, 4, 5])],
    )
}

#[test]
fn disconnected_faces_form_stable_components() {
    let components = surface_components(&two_triangles());
    assert_eq!(components.len(), 2);
    assert_eq!(components[0].faces, vec![0]);
    assert_eq!(components[1].faces, vec![1]);
    assert!(components[0].area > components[1].area);
}

#[test]
fn filtering_keeps_and_compacts_the_largest_component() {
    let result = filter_surface_components(
        &two_triangles(),
        SurfaceComponentFilter {
            minimum_area: 0.0,
            maximum_components: Some(1),
        },
    );
    let Ok(filtered) = result else {
        panic!("component filter should succeed")
    };
    assert_eq!(filtered.faces, vec![SurfaceFace([0, 1, 2])]);
    assert_eq!(filtered.vertices.len(), 3);
}

#[test]
fn a_zero_component_limit_is_rejected() {
    let result = filter_surface_components(
        &two_triangles(),
        SurfaceComponentFilter {
            minimum_area: 0.0,
            maximum_components: Some(0),
        },
    );
    assert_eq!(result, Err(SurfaceComponentError::InvalidFilter));
}
