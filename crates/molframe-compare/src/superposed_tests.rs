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
fn the_fit_maximises_the_tm_sum_rather_than_the_squared_deviation() {
    // Sixty points along a line, then forty more standing 12 Å off it. The
    // least-squares fit splits the difference between the two groups, because
    // it weights a large deviation by its square; the TM kernel weights a point
    // by `1/(1+(d/d0)^2)`, so it prefers a fit that keeps the long group and
    // abandons the short one. The published score is the maximum over
    // superpositions, so it takes the larger of the two.
    let mut reference = Vec::new();
    let mut model = Vec::new();
    for index in 0..60_u16 {
        let x = f32::from(index);
        reference.push([x, 0.0, 0.0]);
        model.push([x, 0.0, 0.0]);
    }
    for index in 0..40_u16 {
        let x = f32::from(index);
        reference.push([x, 5.0, 0.0]);
        model.push([x, 5.0, 12.0]);
    }

    let score = match tm_score(&model, &reference) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    let plain = plain_fit_tm(&model, &reference);
    assert!(
        score > plain,
        "the maximising fit ({score}) must beat the least-squares one ({plain})"
    );
    // Keeping the stationary sixty of a hundred points is worth 0.6 on its own,
    // so a fit that has not found that is still on the compromise solution.
    assert!(score > 0.6, "the stationary group should be kept: {score}");
}

/// The score at the plain least-squares fit, for comparison.
fn plain_fit_tm(model: &[[f32; 3]], reference: &[[f32; 3]]) -> f64 {
    let fit = match molframe_geom::superpose(model, reference) {
        Ok(fit) => fit,
        Err(error) => panic!("valid: {error}"),
    };
    let d0 = 0.5f64;
    let mut total = 0.0f64;
    for (moved, &target) in model
        .iter()
        .map(|&point| fit.transform.apply(point))
        .zip(reference)
    {
        let ratio = molframe_geom::distance(moved, target) / d0;
        total += 1.0 / (1.0 + ratio * ratio);
    }
    let count = match u32::try_from(reference.len()) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    };
    total / f64::from(count)
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

#[test]
fn weighted_rmsd_reports_a_short_weight_vector_as_a_weight_problem() {
    let reference = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let model = [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    assert!(matches!(
        weighted_rmsd(&model, &reference, &[1.0]),
        Err(CompareError::WeightLengthMismatch {
            points: 2,
            weights: 1
        })
    ));
}

#[test]
fn weighted_rmsd_reports_unequal_coordinates_as_a_point_problem() {
    let model = [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    assert!(matches!(
        weighted_rmsd(&model, &[[0.0, 0.0, 0.0]], &[1.0, 1.0]),
        Err(CompareError::LengthMismatch {
            model: 2,
            reference: 1
        })
    ));
}
