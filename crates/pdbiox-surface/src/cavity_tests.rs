use super::{cavities, cavities_with_options};
use crate::{SasaError, SurfaceGridOptions};
use crate::{numeric::f64_to_f32, sampling::fibonacci_sphere};

#[test]
fn a_lone_atom_encloses_nothing() {
    let found = cavities(&[[0.0, 0.0, 0.0]], &[2.0], 0.0, 0.5);
    let Ok(found) = found else { panic!("valid") };
    assert!(found.is_empty());
}

#[test]
fn a_sealed_shell_of_atoms_encloses_one_cavity() {
    // A dense shell on a sphere of radius 5, each atom large enough that its
    // neighbours overlap, seals a hollow interior.
    let positions: Vec<[f32; 3]> = fibonacci_sphere(120)
        .into_iter()
        .map(|direction| {
            [
                f64_to_f32(direction[0] * 5.0),
                f64_to_f32(direction[1] * 5.0),
                f64_to_f32(direction[2] * 5.0),
            ]
        })
        .collect();
    let radii = vec![2.2f32; positions.len()];
    let Ok(found) = cavities(&positions, &radii, 0.0, 0.6) else {
        panic!("valid");
    };
    assert!(!found.is_empty(), "the hollow interior should be a cavity");
    assert!(found[0].volume > 0.0);
    // The largest cavity's representative point sits near the centre.
    let point = found[0].representative;
    let distance =
        (f64::from(point[0]).powi(2) + f64::from(point[1]).powi(2) + f64::from(point[2]).powi(2))
            .sqrt();
    assert!(
        distance < 3.0,
        "cavity point should be interior, was {distance}"
    );
}

#[test]
fn a_non_positive_resolution_is_rejected() {
    assert!(cavities(&[[0.0, 0.0, 0.0]], &[2.0], 0.0, 0.0).is_err());
}

#[test]
fn allocation_ceiling_is_caller_controlled() {
    let error = cavities_with_options(
        &[[0.0, 0.0, 0.0]],
        &[2.0],
        0.0,
        SurfaceGridOptions {
            resolution: 0.5,
            max_cells: 1,
        },
    )
    .expect_err("one cell cannot hold the bounded grid");
    assert!(matches!(error, SasaError::GridTooLarge { .. }));
    assert!(matches!(
        cavities_with_options(
            &[[0.0, 0.0, 0.0]],
            &[2.0],
            0.0,
            SurfaceGridOptions {
                resolution: 0.5,
                max_cells: 0,
            },
        ),
        Err(SasaError::InvalidGridOptions)
    ));
}
