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

/// A deterministic extended C-alpha-like trace with 3.8 A steps.
fn trace(length: usize) -> Vec<[f32; 3]> {
    let mut state = 0x9e37_79b9_u32;
    let mut next = move || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        f32::from(u16::try_from(state >> 16).expect("sixteen bits fit")) / 65_535.0 - 0.5
    };
    let mut points = Vec::with_capacity(length);
    let mut at = [0.0f32; 3];
    for _ in 0..length {
        points.push(at);
        let mut step = [next(), next(), next()];
        let norm = (step[0] * step[0] + step[1] * step[1] + step[2] * step[2]).sqrt();
        for axis in &mut step {
            *axis *= 3.8 / norm;
        }
        at = [at[0] + step[0] + 1.5, at[1] + step[1], at[2] + step[2]];
    }
    points
}

/// Rotates the tail about z through the origin and shifts it.
fn hinge(points: &mut [[f32; 3]], degrees: f32, shift: [f32; 3]) {
    let (sin, cos) = degrees.to_radians().sin_cos();
    for point in points {
        let [x, y, z] = *point;
        *point = [
            cos * x - sin * y + shift[0],
            sin * x + cos * y + shift[1],
            z + shift[2],
        ];
    }
}

#[test]
fn a_hinged_pair_scores_its_best_domain_not_the_global_compromise() {
    // 214 atoms split 55/45, the tail rotated and shifted about 12 A. One domain
    // fits exactly, so the score is at least that domain's share; a single
    // least-squares fit lands far below it.
    let reference = trace(214);
    let mut model = reference.clone();
    hinge(&mut model[118..], 40.0, [12.0, 12.0, 12.0]);
    let score = tm_score(&model, &reference).expect("valid");
    let plain = plain_fit_tm(&model, &reference);
    assert!(score >= 118.0 / 214.0 - 1e-9, "{score}");
    assert!(score > plain + 0.05, "{score} vs plain {plain}");
    assert!(score <= 1.0);
}

#[test]
fn outliers_do_not_hide_a_perfect_core() {
    let reference = trace(200);
    let mut model = reference.clone();
    for point in &mut model[..30] {
        point[2] += 25.0;
    }
    let score = tm_score(&model, &reference).expect("valid");
    assert!(score >= 170.0 / 200.0 - 1e-9, "{score}");
}

#[test]
fn the_tm_score_does_not_depend_on_the_worker_count() {
    let reference = trace(2_500);
    let mut model = reference.clone();
    hinge(&mut model[1_300..], 35.0, [10.0, -8.0, 9.0]);
    let run = |workers: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .expect("pool")
            .install(|| tm_score(&model, &reference).expect("valid"))
    };
    let baseline = run(1);
    for workers in [2, 5] {
        assert_eq!(baseline.to_bits(), run(workers).to_bits(), "{workers}");
    }
}

#[test]
fn twenty_thousand_atoms_score_in_bounded_time() {
    let reference = trace(20_000);
    let mut model = reference.clone();
    hinge(&mut model[11_000..], 40.0, [30.0, 20.0, 25.0]);
    let started = std::time::Instant::now();
    let score = tm_score(&model, &reference).expect("valid");
    assert!(started.elapsed().as_secs() < 60, "{:?}", started.elapsed());
    assert!((0.55..=1.0).contains(&score), "{score}");
}
