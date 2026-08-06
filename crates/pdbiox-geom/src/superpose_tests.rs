#![allow(
    clippy::float_cmp,
    reason = "exact identities on exactly-representable inputs"
)]

use super::*;

const TRIANGLE: [[f32; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

#[test]
fn a_set_fits_itself_exactly_and_the_transform_changes_nothing() {
    let Ok(fit) = superpose(&TRIANGLE, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(fit.rmsd < 1e-6, "got {}", fit.rmsd);
    for point in TRIANGLE {
        let moved = fit.transform.apply(point);
        for axis in 0..3 {
            assert!((moved[axis] - point[axis]).abs() < 1e-5);
        }
    }
}

#[test]
fn a_translated_copy_fits_exactly() {
    let moved: Vec<[f32; 3]> = TRIANGLE
        .iter()
        .map(|p| [p[0] + 10.0, p[1] - 5.0, p[2] + 2.0])
        .collect();
    let Ok(fit) = superpose(&moved, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(fit.rmsd < 1e-5, "got {}", fit.rmsd);
}

#[test]
fn a_rotated_copy_fits_exactly_and_the_fit_is_a_rotation_not_a_reflection() {
    // A quarter turn about z.
    let turned: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [-p[1], p[0], p[2]]).collect();
    let Ok(fit) = superpose(&turned, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(fit.rmsd < 1e-5, "got {}", fit.rmsd);
    assert!(
        (fit.transform.determinant() - 1.0).abs() < 1e-9,
        "a fit must never mirror the structure"
    );
}

#[test]
fn a_mirrored_copy_is_not_fitted_by_reflecting_it() {
    let mirrored: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [-p[0], p[1], p[2]]).collect();
    let Ok(fit) = superpose(&mirrored, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(
        (fit.transform.determinant() - 1.0).abs() < 1e-9,
        "the transform must stay a rotation"
    );
    assert!(
        fit.rmsd > 0.1,
        "a mirror image genuinely does not fit, got {}",
        fit.rmsd
    );
}

#[test]
fn deviation_is_symmetric_and_zero_against_itself() {
    assert_eq!(rmsd(&TRIANGLE, &TRIANGLE).ok(), Some(0.0));
    let other: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [p[0] + 1.0, p[1], p[2]]).collect();
    let (Ok(forward), Ok(backward)) = (rmsd(&TRIANGLE, &other), rmsd(&other, &TRIANGLE)) else {
        panic!("expected both")
    };
    assert!((forward - backward).abs() < 1e-12);
    assert!((forward - 1.0).abs() < 1e-6);
}

#[test]
fn sets_of_different_sizes_are_refused_rather_than_truncated() {
    assert_eq!(
        rmsd(&TRIANGLE, &TRIANGLE[..2]),
        Err(SuperposeError::LengthMismatch)
    );
    assert_eq!(
        superpose(&TRIANGLE, &TRIANGLE[..2]).err(),
        Some(SuperposeError::LengthMismatch)
    );
}

#[test]
fn too_few_points_to_fix_a_rotation_are_refused() {
    let pair = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    assert_eq!(
        superpose(&pair, &pair).err(),
        Some(SuperposeError::TooFewPoints)
    );
}

#[test]
fn the_same_pair_of_sets_always_fits_the_same_way() {
    let moved: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [p[1], p[2], p[0]]).collect();
    let (Ok(first), Ok(second)) = (superpose(&moved, &TRIANGLE), superpose(&moved, &TRIANGLE))
    else {
        panic!("expected both")
    };
    assert_eq!(first.transform, second.transform);
    assert_eq!(first.rmsd, second.rmsd);
}
