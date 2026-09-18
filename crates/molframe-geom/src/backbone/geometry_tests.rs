use super::{backbone_frames, helix_geometry};
use std::f64::consts::FRAC_PI_2;

#[test]
fn missing_points_invalidate_only_their_local_frames() {
    let points = [
        Some([0.0, 0.0, 0.0]),
        Some([1.0, 0.0, 0.0]),
        None,
        Some([2.0, 1.0, 0.0]),
        Some([3.0, 1.0, 0.0]),
    ];
    assert!(backbone_frames(&points).iter().all(Option::is_none));
}

#[test]
fn a_right_angle_has_the_expected_discrete_curvature() {
    let frames = backbone_frames(&[
        Some([0.0, 0.0, 0.0]),
        Some([1.0, 0.0, 0.0]),
        Some([1.0, 1.0, 0.0]),
    ]);
    assert!(frames[1].is_some_and(|frame| (frame.curvature - FRAC_PI_2).abs() < 1e-12));
}

#[test]
fn a_regular_helix_reports_positive_rise_and_quarter_turn_twist() {
    let points = [
        Some([1.0, 0.0, 0.0]),
        Some([0.0, 1.0, 1.5]),
        Some([-1.0, 0.0, 3.0]),
        Some([0.0, -1.0, 4.5]),
        Some([1.0, 0.0, 6.0]),
        Some([0.0, 1.0, 7.5]),
        Some([-1.0, 0.0, 9.0]),
        Some([0.0, -1.0, 10.5]),
    ];
    let geometry = helix_geometry(&points);
    assert!(geometry.is_ok_and(|geometry| {
        geometry
            .is_some_and(|value| value.rise > 1.4 && (value.twist.abs() - FRAC_PI_2).abs() < 0.2)
    }));
}
