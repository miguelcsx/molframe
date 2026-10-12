//! Deviation and flat-lane agreement.

use super::{rmsd, rmsd_flat};
use crate::SuperposeError;
use crate::superpose::superpose;

const TRIANGLE: [[f32; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

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
fn flat_coordinates_use_the_same_kernel_without_truncation() {
    let flat = TRIANGLE.into_iter().flatten().collect::<Vec<_>>();
    assert_eq!(rmsd_flat(&flat, &flat), Ok(0.0));
    assert_eq!(
        rmsd_flat(&flat[..flat.len() - 1], &flat[..flat.len() - 1]),
        Err(SuperposeError::LengthMismatch)
    );
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
