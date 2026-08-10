use super::{gdt_ts, gdt_with_cutoffs, tm_score, weighted_rmsd};
use crate::failure::CompareError;

const TETRAHEDRON: [[f32; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn translated(offset: f32) -> [[f32; 3]; 4] {
    let mut moved = TETRAHEDRON;
    for point in &mut moved {
        point[0] += offset;
        point[1] += offset;
        point[2] += offset;
    }
    moved
}

#[test]
fn an_identical_structure_scores_one_on_both_measures() {
    let tm = match tm_score(&TETRAHEDRON, &TETRAHEDRON) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    let gdt = match gdt_ts(&TETRAHEDRON, &TETRAHEDRON) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!((tm - 1.0).abs() < 1e-6, "tm {tm}");
    assert!((gdt - 1.0).abs() < 1e-6, "gdt {gdt}");
}

#[test]
fn a_translated_copy_is_fitted_back_to_a_perfect_score() {
    let model = translated(10.0);
    let tm = match tm_score(&model, &TETRAHEDRON) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!((tm - 1.0).abs() < 1e-5, "translation should fit out: {tm}");
}

#[test]
fn unequal_lengths_are_rejected() {
    let short = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(matches!(
        tm_score(&TETRAHEDRON, &short),
        Err(CompareError::LengthMismatch { .. })
    ));
}

#[test]
fn too_few_points_cannot_be_superposed() {
    let two = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    assert!(matches!(
        tm_score(&two, &two),
        Err(CompareError::Superpose(_))
    ));
}

#[test]
fn general_gdt_kernel_requires_explicit_ordered_cutoffs() {
    let score = gdt_with_cutoffs(&TETRAHEDRON, &TETRAHEDRON, &[0.25, 0.75])
        .unwrap_or_else(|error| panic!("explicit GDT failed: {error}"));
    assert!((score - 1.0).abs() < 1.0e-12);
    assert!(matches!(
        gdt_with_cutoffs(&TETRAHEDRON, &TETRAHEDRON, &[2.0, 1.0]),
        Err(CompareError::InvalidDistanceCutoff)
    ));
}

#[test]
fn weighted_rmsd_uses_only_explicit_nonnegative_weights() {
    let reference = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let model = [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let Ok(value) = weighted_rmsd(&model, &reference, &[3.0, 1.0]) else {
        panic!("valid weighted correspondence");
    };
    assert!((value - 3.0_f64.sqrt()).abs() < 1.0e-12);
    assert!(matches!(
        weighted_rmsd(&model, &reference, &[0.0, 0.0]),
        Err(CompareError::InvalidScoreInput)
    ));
}
