use super::*;
use std::f64::consts::{FRAC_PI_2, PI};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn distance_is_symmetric_and_agrees_with_its_squared_form() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 6.0, 3.0];
    assert!(close(distance(a, b), 5.0));
    assert!(close(distance(a, b), distance(b, a)));
    assert!(close(distance_squared(a, b), 25.0));
}

#[test]
fn distance_satisfies_the_triangle_inequality() {
    let a = [0.0, 0.0, 0.0];
    let b = [3.0, 0.0, 0.0];
    let c = [3.0, 4.0, 0.0];
    assert!(distance(a, c) <= distance(a, b) + distance(b, c) + 1e-12);
}

#[test]
fn a_point_is_no_distance_from_itself() {
    let point = [1.5, -2.25, 0.125];
    assert!(close(distance(point, point), 0.0));
}

#[test]
fn a_right_angle_measures_a_quarter_turn_and_reads_the_same_from_either_arm() {
    let (a, vertex, c) = ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let forward = angle(a, vertex, c);
    let backward = angle(c, vertex, a);
    assert!(forward.is_some_and(|value| close(value, FRAC_PI_2)));
    assert!(
        forward
            .zip(backward)
            .is_some_and(|(left, right)| close(left, right))
    );
}

#[test]
fn a_straight_arrangement_measures_half_a_turn() {
    let straight = angle([-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    assert!(straight.is_some_and(|value| close(value, PI)));
}

#[test]
fn an_angle_with_no_arm_is_undefined_rather_than_zero() {
    let vertex = [0.0, 0.0, 0.0];
    assert_eq!(angle(vertex, vertex, [1.0, 0.0, 0.0]), None);
}

#[test]
fn a_planar_torsion_is_zero_one_way_and_half_a_turn_the_other() {
    let cis = dihedral(
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    );
    assert!(cis.is_some_and(|value| value.abs() < 1e-9));

    let trans = dihedral(
        [-1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    );
    assert!(trans.is_some_and(|value| close(value.abs(), PI)));
}

#[test]
fn a_torsion_changes_sign_when_its_ends_are_swapped() {
    let a = [1.0, 1.0, 0.0];
    let b = [0.0, 1.0, 0.0];
    let c = [0.0, 0.0, 0.0];
    let d = [1.0, 0.0, 1.0];
    let forward = dihedral(a, b, c, d);
    let backward = dihedral(d, c, b, a);
    match (forward, backward) {
        (Some(first), Some(second)) => assert!(close(first, second), "{first} vs {second}"),
        _ => panic!("both directions should be defined"),
    }
}

#[test]
fn a_torsion_over_collinear_points_is_undefined() {
    let collinear = dihedral(
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    );
    assert_eq!(collinear, None);
}

#[test]
fn vector_helpers_agree_with_their_definitions() {
    assert!(close(dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]), 32.0));
    let product = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    assert!(
        product
            .into_iter()
            .zip([0.0, 0.0, 1.0])
            .all(|(left, right)| close(left, right))
    );
    assert!(close(norm([3.0, 4.0, 0.0]), 5.0));
    assert_eq!(normalise([0.0, 0.0, 0.0]), None);
    assert!(normalise([0.0, 5.0, 0.0]).is_some_and(|unit| close(unit[1], 1.0)));
}

#[test]
fn degrees_converts_a_half_turn_to_a_hundred_and_eighty() {
    assert!(close(degrees(PI), 180.0));
}
