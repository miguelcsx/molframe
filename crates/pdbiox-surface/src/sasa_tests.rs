use super::{SasaError, shrake_rupley, surface_points};
use core::f64::consts::PI;

fn sphere_area(radius: f64) -> f64 {
    4.0 * PI * radius * radius
}

// The kernel widens each radius with the probe in f64 after an f32 input, so a
// reference area must expand the same way rather than assume exact decimals.
fn expanded(radius: f32, probe: f32) -> f64 {
    f64::from(radius) + f64::from(probe)
}

#[test]
fn a_lone_atom_exposes_its_whole_expanded_sphere() {
    let Ok(areas) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.5], 1.4, 400) else {
        panic!("a single atom is valid");
    };
    assert!((areas[0] - sphere_area(expanded(1.5, 1.4))).abs() < 1e-6);
}

#[test]
fn a_small_atom_buried_inside_a_larger_one_has_no_accessible_area() {
    // Concentric: the small sphere lies wholly within the large one.
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let Ok(areas) = shrake_rupley(&positions, &[1.0, 3.0], 0.0, 400) else {
        panic!("valid");
    };
    assert!(areas[0].abs() < 1e-9, "buried atom kept area {}", areas[0]);
    assert!((areas[1] - sphere_area(3.0)).abs() < 1e-6);
}

#[test]
fn atoms_far_apart_do_not_shade_each_other() {
    let positions = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
    let Ok(areas) = shrake_rupley(&positions, &[1.5, 1.5], 1.4, 400) else {
        panic!("valid");
    };
    let full = sphere_area(expanded(1.5, 1.4));
    assert!((areas[0] - full).abs() < 1e-6);
    assert!((areas[1] - full).abs() < 1e-6);
}

#[test]
fn touching_atoms_lose_area_and_share_the_loss() {
    // Two equal spheres whose expanded radii overlap occlude each other equally.
    // The physical areas are identical; a finite, non-mirror-symmetric point
    // lattice only reproduces that to within its sampling, so the two match
    // closely rather than exactly, and both sit below a full sphere.
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let Ok(areas) = shrake_rupley(&positions, &[1.5, 1.5], 0.0, 2000) else {
        panic!("valid");
    };
    let full = sphere_area(1.5);
    assert!(
        (areas[0] - areas[1]).abs() < 0.01 * full,
        "occlusion should be near-symmetric: {} vs {}",
        areas[0],
        areas[1]
    );
    assert!(areas[0] < full && areas[0] > 0.0, "expected partial area");
    assert!(areas[1] < full && areas[1] > 0.0, "expected partial area");
}

#[test]
fn a_mismatched_radius_count_is_rejected() {
    let Err(error) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.0, 2.0], 1.4, 100) else {
        panic!("mismatch should be rejected");
    };
    assert!(matches!(error, SasaError::LengthMismatch { .. }));
}

#[test]
fn a_negative_probe_is_rejected() {
    let Err(error) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.0], -1.0, 100) else {
        panic!("negative probe should be rejected");
    };
    assert!(matches!(error, SasaError::InvalidProbe));
}

#[test]
fn zero_test_points_is_rejected() {
    let Err(error) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.0], 1.4, 0) else {
        panic!("zero points should be rejected");
    };
    assert!(matches!(error, SasaError::NoPoints));
}

#[test]
fn an_empty_set_has_no_areas() {
    let Ok(areas) = shrake_rupley(&[], &[], 1.4, 100) else {
        panic!("an empty set is valid");
    };
    assert!(areas.is_empty());
}

#[test]
fn a_lone_atom_yields_one_surface_point_per_sample() {
    let Ok(points) = surface_points(&[[0.0, 0.0, 0.0]], &[1.5], 1.4, 200) else {
        panic!("valid");
    };
    assert_eq!(points.len(), 200);
    for point in &points {
        assert_eq!(point.atom, 0);
        // Every point sits at the expanded radius from the centre.
        let radius = expanded(1.5, 1.4);
        let distance = (f64::from(point.position[0]).powi(2)
            + f64::from(point.position[1]).powi(2)
            + f64::from(point.position[2]).powi(2))
        .sqrt();
        assert!((distance - radius).abs() < 1e-4, "distance {distance}");
    }
}

#[test]
fn a_buried_atom_contributes_no_surface_points() {
    // A small atom concentric inside a large one is fully buried.
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let Ok(points) = surface_points(&positions, &[1.0, 3.0], 0.0, 400) else {
        panic!("valid");
    };
    assert!(
        points.iter().all(|point| point.atom != 0),
        "atom 0 is buried"
    );
}
