use super::lee_richards;
use crate::accessible_area::shrake_rupley;
use core::f64::consts::PI;
use molframe_core::ExecutionContext;

fn expanded(radius: f32, probe: f32) -> f64 {
    f64::from(radius) + f64::from(probe)
}

#[test]
fn a_lone_atom_recovers_its_sphere_area_exactly() {
    let Ok(areas) = lee_richards(
        &[[0.0, 0.0, 0.0]],
        &[1.5],
        1.4,
        300,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    let radius = expanded(1.5, 1.4);
    assert!(
        (areas[0] - 4.0 * PI * radius * radius).abs() < 1e-6,
        "area {}",
        areas[0]
    );
}

#[test]
fn far_apart_atoms_are_each_fully_exposed() {
    let positions = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
    let Ok(areas) = lee_richards(
        &positions,
        &[1.5, 1.5],
        1.4,
        300,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    let full = 4.0 * PI * expanded(1.5, 1.4).powi(2);
    assert!((areas[0] - full).abs() < 1e-6);
    assert!((areas[1] - full).abs() < 1e-6);
}

#[test]
fn a_buried_atom_has_almost_no_accessible_area() {
    // A small atom concentric within a much larger one is essentially buried.
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let Ok(areas) = lee_richards(
        &positions,
        &[1.0, 3.0],
        0.0,
        400,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    let small_full = 4.0 * PI;
    assert!(
        areas[0] < 0.02 * small_full,
        "buried atom kept {}",
        areas[0]
    );
}

#[test]
fn it_agrees_with_shrake_rupley_on_a_touching_pair() {
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let Ok(lr) = lee_richards(
        &positions,
        &[1.5, 1.5],
        1.4,
        400,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    let Ok(sr) = shrake_rupley(
        &positions,
        &[1.5, 1.5],
        1.4,
        4000,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    let lr_total: f64 = lr.iter().sum();
    let sr_total: f64 = sr.iter().sum();
    assert!(
        (lr_total - sr_total).abs() < 0.05 * sr_total,
        "lee-richards {lr_total} vs shrake-rupley {sr_total}"
    );
}

#[test]
fn zero_slices_is_rejected() {
    assert!(
        lee_richards(
            &[[0.0, 0.0, 0.0]],
            &[1.5],
            1.4,
            0,
            &ExecutionContext::default(),
        )
        .is_err()
    );
}
