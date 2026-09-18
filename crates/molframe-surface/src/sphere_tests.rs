use super::fibonacci_sphere;
use crate::numeric::usize_to_f64;

#[test]
fn every_direction_is_a_unit_vector() {
    for direction in fibonacci_sphere(200) {
        let length = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        assert!((length - 1.0).abs() < 1e-9, "off the unit sphere: {length}");
    }
}

#[test]
fn the_lattice_is_balanced_so_the_centroid_sits_near_the_origin() {
    let points = fibonacci_sphere(1000);
    let mut sum = [0.0f64; 3];
    for point in &points {
        sum[0] += point[0];
        sum[1] += point[1];
        sum[2] += point[2];
    }
    let mean =
        (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt() / usize_to_f64(points.len());
    assert!(mean < 1e-2, "lattice is lopsided: {mean}");
}

#[test]
fn the_same_count_returns_identical_directions() {
    assert_eq!(fibonacci_sphere(64), fibonacci_sphere(64));
}

#[test]
fn a_request_for_no_points_is_empty() {
    assert!(fibonacci_sphere(0).is_empty());
}
