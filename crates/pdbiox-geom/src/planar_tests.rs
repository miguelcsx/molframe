use super::{best_fit_plane, plane_deviation};

#[test]
fn fewer_than_three_points_have_no_plane() {
    assert!(
        plane_deviation(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            .is_ok_and(|deviation| deviation.is_none())
    );
}

#[test]
fn coplanar_points_deviate_by_nothing() {
    let square = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let Ok(Some(deviation)) = plane_deviation(&square) else {
        panic!("four points define a plane");
    };
    assert!(deviation < 1e-6, "deviation {deviation}");
}

#[test]
fn a_tilted_plane_is_still_flat() {
    // Points on the plane x = y still lie in a plane, just not an axis-aligned one.
    let tilted = [
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 3.0],
        [0.0, 0.0, 2.0],
    ];
    let Ok(Some(deviation)) = plane_deviation(&tilted) else {
        panic!("valid");
    };
    assert!(deviation < 1e-6, "deviation {deviation}");
}

#[test]
fn a_puckered_set_departs_from_any_plane() {
    let tetrahedron = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let Ok(Some(deviation)) = plane_deviation(&tetrahedron) else {
        panic!("valid");
    };
    assert!(deviation > 0.1, "a tetrahedron is not planar: {deviation}");
}

#[test]
fn the_normal_of_an_xy_plane_points_along_z() {
    let square = [
        [0.0, 0.0, 5.0],
        [1.0, 0.0, 5.0],
        [1.0, 1.0, 5.0],
        [0.0, 1.0, 5.0],
    ];
    let Ok(Some(plane)) = best_fit_plane(&square) else {
        panic!("four points define a plane");
    };
    // The centroid sits in the plane, and the normal is ±z.
    assert!((plane.centre[2] - 5.0).abs() < 1e-9);
    assert!(plane.normal[0].abs() < 1e-6 && plane.normal[1].abs() < 1e-6);
    assert!((plane.normal[2].abs() - 1.0).abs() < 1e-6);
}
