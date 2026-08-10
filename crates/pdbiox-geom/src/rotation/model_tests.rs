use super::{Rotation3, RotationError, RotationOptions, rotation_mean};
use std::f64::consts::{FRAC_PI_2, PI};

fn around_z(angle: f64) -> Rotation3 {
    match Rotation3::exp([0.0, 0.0, angle]) {
        Ok(rotation) => rotation,
        Err(error) => panic!("valid rotation rejected: {error:?}"),
    }
}

#[test]
fn reflections_and_shears_are_refused() {
    let reflection = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        Rotation3::from_matrix(reflection, 1e-12),
        Err(RotationError::Reflection)
    );
    let shear = [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        Rotation3::from_matrix(shear, 1e-12),
        Err(RotationError::NotOrthonormal)
    );
}

#[test]
fn exponential_and_logarithm_round_trip() {
    let rotation = around_z(0.7);
    let rebuilt = Rotation3::exp(rotation.log());
    assert!(rebuilt.is_ok_and(|rebuilt| rebuilt.distance(rotation) < 1e-10));
}

#[test]
fn geodesic_distance_is_invariant_under_common_rotation() {
    let first = around_z(0.2);
    let second = around_z(1.1);
    let common = match Rotation3::exp([0.4, 0.2, -0.1]) {
        Ok(value) => value,
        Err(error) => panic!("valid rotation rejected: {error:?}"),
    };
    assert!(
        (first.distance(second) - first.then(common).distance(second.then(common))).abs() < 1e-10
    );
}

#[test]
fn interpolation_follows_constant_speed() {
    let start = Rotation3::IDENTITY;
    let end = around_z(PI);
    let middle = start.interpolate(end, 0.5);
    assert!(middle.is_ok_and(|middle| (start.distance(middle) - FRAC_PI_2).abs() < 1e-8));
}

#[test]
fn intrinsic_mean_sits_between_nearby_rotations() {
    let mean = rotation_mean(&[around_z(-0.2), around_z(0.2)], 1e-12, 32);
    assert!(mean.is_ok_and(|mean| mean.distance(Rotation3::IDENTITY) < 1e-10));
}

#[test]
fn inconsistent_rotation_controls_are_refused() {
    let options = RotationOptions {
        matrix_tolerance: 1e-12,
        small_angle_tolerance: 1e-2,
        half_turn_tolerance: 1e-8,
    };
    assert_eq!(
        Rotation3::exp_with_options([0.0, 0.0, 1e-3], options),
        Err(RotationError::InvalidOptions)
    );
}
