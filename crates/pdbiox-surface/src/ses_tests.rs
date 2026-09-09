use super::{solvent_excluded_surface, solvent_excluded_surface_with_options};
use crate::{SasaError, SurfaceGridOptions};

#[test]
fn a_lone_atom_converges_on_its_vdw_sphere_not_the_grown_sphere() {
    let surface = solvent_excluded_surface(&[[0.0, 0.0, 0.0]], &[2.0], 1.4, 0.1)
        .unwrap_or_else(|error| panic!("SES failed: {error}"));
    let expected = 4.0 * core::f64::consts::PI * 4.0;
    assert!(
        (surface.area - expected).abs() / expected < 0.2,
        "area {}",
        surface.area
    );
    assert!(!surface.triangles.is_empty());
    assert!(surface.triangles.iter().all(|triangle| triangle.area > 0.0));
}

#[test]
fn empty_input_has_an_empty_surface() {
    let surface = solvent_excluded_surface(&[], &[], 1.4, 0.25)
        .unwrap_or_else(|error| panic!("SES failed: {error}"));
    assert!(surface.triangles.is_empty());
    assert!(surface.area.abs() < f64::EPSILON);
}

#[test]
fn a_touching_pair_forms_one_surface_with_a_reentrant_neck() {
    let surface =
        solvent_excluded_surface(&[[0.0, 0.0, 0.0], [3.2, 0.0, 0.0]], &[2.0, 2.0], 1.4, 0.2)
            .unwrap_or_else(|error| panic!("SES failed: {error}"));
    let two_isolated_spheres = 2.0 * 4.0 * core::f64::consts::PI * 4.0;
    assert!(surface.area > 0.0 && surface.area < two_isolated_spheres);
    assert!(surface.triangles.iter().any(|triangle| {
        triangle
            .vertices
            .iter()
            .any(|vertex| vertex[0] > 1.0 && vertex[0] < 2.2)
    }));
}

#[test]
fn invalid_resolution_is_rejected() {
    assert!(solvent_excluded_surface(&[[0.0; 3]], &[1.0], 1.4, 0.0).is_err());
}

#[test]
fn allocation_ceiling_is_caller_controlled() {
    let result = solvent_excluded_surface_with_options(
        &[[0.0; 3]],
        &[1.0],
        1.4,
        SurfaceGridOptions {
            resolution: 0.5,
            max_cells: 1,
            max_workspace_bytes: SurfaceGridOptions::STANDARD_WORKSPACE_BYTES,
        },
    );

    assert!(matches!(result, Err(SasaError::GridTooLarge { .. })));
}
