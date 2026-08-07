use super::place_atom;
use crate::Dihedron;

#[test]
fn measured_geometry_places_the_fourth_point_back_where_it_started() {
    let points = [
        [0.2, -0.1, 0.3],
        [1.1, 0.0, 0.2],
        [1.4, 1.2, -0.2],
        [2.0, 1.5, 0.9],
    ];
    let Some(internal) = Dihedron::from_points(points[0], points[1], points[2], points[3]) else {
        panic!("geometry degenerate")
    };
    let Some(rebuilt) = place_atom(
        points[0],
        points[1],
        points[2],
        internal.length,
        internal.angle,
        internal.torsion,
    ) else {
        panic!("placement failed")
    };
    assert!(
        rebuilt
            .iter()
            .zip(points[3])
            .all(|(a, b)| (*a - b).abs() < 1e-5)
    );
}

#[test]
fn collinear_references_are_explicitly_undefined() {
    assert!(
        place_atom(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            1.0,
            1.0,
            0.0,
        )
        .is_none()
    );
}
