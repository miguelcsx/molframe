use super::{PeriodicAngle, PeriodicError, TorusMetric, circular_summary, torus_summary};
use std::f64::consts::{PI, TAU};

fn angle(value: f64) -> PeriodicAngle {
    match PeriodicAngle::from_radians(value) {
        Ok(angle) => angle,
        Err(error) => panic!("finite angle rejected: {error:?}"),
    }
}

#[test]
fn whole_turns_are_the_same_angle() {
    assert!((angle(0.3).signed_delta(angle(0.3 + 4.0 * TAU))).abs() < 1e-12);
}

#[test]
fn the_half_turn_tie_has_one_canonical_sign() {
    assert!((angle(-PI).radians() - PI).abs() < f64::EPSILON);
    assert!((angle(PI).radians() - PI).abs() < f64::EPSILON);
}

#[test]
fn circular_distance_crosses_the_branch_cut() {
    let left = angle(179_f64.to_radians());
    let right = angle(-179_f64.to_radians());
    assert!((left.distance(right).to_degrees() - 2.0).abs() < 1e-10);
}

#[test]
fn antipodal_samples_have_no_defensible_mean() {
    let summary = circular_summary(&[angle(0.0), angle(PI)], 1e-12);
    let Ok(summary) = summary else {
        panic!("valid sample rejected");
    };
    assert!(summary.mean.is_none());
    assert!(summary.resultant < 1e-12);
}

#[test]
fn torus_distance_is_wrapped_weighted_and_symmetric() {
    let metric = TorusMetric::new(vec![1.0, 4.0].into_boxed_slice());
    let Ok(metric) = metric else {
        panic!("valid metric rejected");
    };
    let left = [angle(179_f64.to_radians()), angle(0.0)];
    let right = [angle(-179_f64.to_radians()), angle(0.5)];
    let forward = metric.distance(&left, &right);
    let reverse = metric.distance(&right, &left);
    assert_eq!(forward, reverse);
    assert!(forward.is_ok_and(|value| value < 1.01));
}

#[test]
fn a_metric_refuses_invalid_weights_and_dimensions() {
    assert_eq!(
        TorusMetric::new(vec![0.0].into_boxed_slice()),
        Err(PeriodicError::InvalidWeight)
    );
    let metric = TorusMetric::uniform(2);
    let Ok(metric) = metric else {
        panic!("uniform metric rejected");
    };
    assert_eq!(
        metric.distance(&[angle(0.0)], &[angle(0.0)]),
        Err(PeriodicError::DimensionMismatch)
    );
}

#[test]
fn torus_summary_preserves_each_component() {
    let points = vec![
        vec![angle(0.1), angle(-0.2)].into_boxed_slice(),
        vec![angle(0.3), angle(0.2)].into_boxed_slice(),
    ];
    let summary = torus_summary(&points, 1e-12);
    assert!(summary.is_ok_and(
        |summary| summary.len() == 2 && summary.iter().all(|value| value.mean.is_some())
    ));
}
