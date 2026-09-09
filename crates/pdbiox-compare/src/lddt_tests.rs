use super::{EmptyLddtPolicy, LddtOptions, lddt, lddt_with_options};
use crate::failure::CompareError;
use pdbiox_core::ExecutionContext;

fn execution(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("valid execution context")
}

#[test]
fn an_identical_model_reproduces_every_distance() {
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
    let score = match lddt(&reference, &reference, 15.0, &execution(1)) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!((score - 1.0).abs() < 1e-12);
}

#[test]
fn a_distance_broken_past_every_tolerance_scores_zero() {
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [13.0, 0.0, 0.0]];
    let score = match lddt(&model, &reference, 15.0, &execution(1)) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!(score.abs() < 1e-12, "score was {score}");
}

#[test]
fn a_distance_off_by_one_and_a_half_keeps_half_the_tolerances() {
    // Reference distance 3, model distance 4.5: an error of 1.5 satisfies the
    // 2 Å and 4 Å tolerances but not 0.5 Å or 1 Å, so the pair scores 0.5.
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [4.5, 0.0, 0.0]];
    let score = match lddt(&model, &reference, 15.0, &execution(1)) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!((score - 0.5).abs() < 1e-9, "score was {score}");
}

#[test]
fn distances_beyond_the_inclusion_radius_do_not_count() {
    // The only pair is 3 Å apart; a 2 Å radius excludes it, leaving nothing to
    // score, which is treated as a perfect (vacuous) result.
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0]];
    let score = match lddt(&model, &reference, 2.0, &execution(1)) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!((score - 1.0).abs() < 1e-12);
}

#[test]
fn unequal_lengths_are_rejected() {
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0]];
    assert!(matches!(
        lddt(&model, &reference, 15.0, &execution(1)),
        Err(CompareError::LengthMismatch { .. })
    ));
}

#[test]
fn general_kernel_uses_explicit_tolerances_and_empty_policy() {
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [4.5, 0.0, 0.0]];
    let score = lddt_with_options(
        &model,
        &reference,
        &LddtOptions {
            inclusion_radius: 10.0,
            minimum_reference_distance: 0.0,
            tolerances: Box::new([1.0, 2.0]),
            empty_policy: EmptyLddtPolicy::Error,
        },
        &execution(1),
    )
    .unwrap_or_else(|error| panic!("explicit lDDT failed: {error}"));
    assert!((score - 0.5).abs() < 1.0e-12);

    assert!(matches!(
        lddt_with_options(
            &model,
            &reference,
            &LddtOptions {
                inclusion_radius: 2.0,
                minimum_reference_distance: 0.0,
                tolerances: Box::new([1.0]),
                empty_policy: EmptyLddtPolicy::Error,
            },
            &execution(1)
        ),
        Err(CompareError::NoComparablePairs)
    ));
}

#[test]
fn indexed_scoring_matches_a_direct_pair_enumeration() {
    let reference: Vec<[f32; 3]> = (0_u16..257)
        .map(|index| {
            let value = f32::from(index);
            [
                value.mul_add(1.37, (value * 0.17).sin()),
                (value * 0.31).cos() * 7.0,
                (value * 0.11).sin() * 3.0,
            ]
        })
        .collect();
    let model: Vec<[f32; 3]> = reference
        .iter()
        .zip(0_u16..)
        .map(|(&[x, y, z], index)| {
            let offset = f32::from(index % 11) * 0.013;
            [x + offset, y - offset * 0.5, z + offset * 0.25]
        })
        .collect();
    let options = LddtOptions {
        inclusion_radius: 12.345_678_9,
        minimum_reference_distance: 0.75,
        tolerances: Box::new([0.05, 0.2, 0.7, 1.5]),
        empty_policy: EmptyLddtPolicy::Error,
    };
    let indexed = lddt_with_options(&model, &reference, &options, &execution(1))
        .unwrap_or_else(|error| panic!("indexed lDDT failed: {error}"));
    let direct = direct_score(&model, &reference, &options);
    assert!((indexed - direct).abs() < 1.0e-15, "{indexed} != {direct}");
}

#[test]
fn a_long_sparse_chain_uses_the_local_domain() {
    let reference: Vec<[f32; 3]> = (0_u16..10_000)
        .map(|index| [f32::from(index) * 3.8, 0.0, 0.0])
        .collect();
    let score = lddt(&reference, &reference, 15.0, &execution(1))
        .unwrap_or_else(|error| panic!("large local lDDT failed: {error}"));
    assert!((score - 1.0).abs() < 1.0e-12);
}

#[test]
fn candidate_rounding_does_not_drop_a_pair_on_the_exact_radius() {
    let reference = [
        [0.661_256_2, 0.268_142_76, 0.099_061_5],
        [0.102_757_365, 0.806_098_4, 0.605_264_9],
    ];
    let radius = test_squared_distance(reference[0], reference[1]).sqrt();
    let options = LddtOptions {
        inclusion_radius: radius,
        minimum_reference_distance: 0.0,
        tolerances: Box::new([0.5, 1.0, 2.0, 4.0]),
        empty_policy: EmptyLddtPolicy::Error,
    };
    assert!(matches!(
        lddt_with_options(&reference, &reference, &options, &execution(1)),
        Ok(score) if (score - 1.0).abs() < f64::EPSILON
    ));
}

fn direct_score(model: &[[f32; 3]], reference: &[[f32; 3]], options: &LddtOptions) -> f64 {
    let minimum_squared = options.minimum_reference_distance.powi(2);
    let inclusion_squared = options.inclusion_radius.powi(2);
    let mut considered = 0_u64;
    let mut preserved = vec![0_u64; options.tolerances.len()];
    for left in 0..reference.len() {
        for right in (left + 1)..reference.len() {
            let squared = test_squared_distance(reference[left], reference[right]);
            if squared <= minimum_squared || squared > inclusion_squared {
                continue;
            }
            considered += 1;
            let model_distance = test_squared_distance(model[left], model[right]).sqrt();
            let error = (model_distance - squared.sqrt()).abs();
            for (count, tolerance) in preserved.iter_mut().zip(&options.tolerances) {
                if error < *tolerance {
                    *count += 1;
                }
            }
        }
    }
    let numerator: u64 = preserved.into_iter().sum();
    crate::numeric::u64_to_f64(numerator)
        / (crate::numeric::u64_to_f64(considered)
            * crate::numeric::usize_to_f64(options.tolerances.len()))
}

fn test_squared_distance(left: [f32; 3], right: [f32; 3]) -> f64 {
    let dx = f64::from(left[0]) - f64::from(right[0]);
    let dy = f64::from(left[1]) - f64::from(right[1]);
    let dz = f64::from(left[2]) - f64::from(right[2]);
    dx.mul_add(dx, dy.mul_add(dy, dz * dz))
}

#[test]
fn the_score_is_identical_at_every_worker_count() {
    // A lattice large enough to span several cell blocks, so the parallel path
    // is genuinely exercised rather than collapsing to one accumulator.
    let reference: Vec<[f32; 3]> = (0..8i16)
        .flat_map(|x| {
            (0..8i16).flat_map(move |y| {
                (0..8i16).map(move |z| [f32::from(x), f32::from(y), f32::from(z)])
            })
        })
        .collect();
    let model: Vec<[f32; 3]> = reference
        .iter()
        .enumerate()
        .map(|(index, point)| {
            let nudge = if index % 3 == 0 { 0.15 } else { -0.08 };
            [point[0] + nudge, point[1], point[2] - nudge]
        })
        .collect();

    let options = LddtOptions::standard(4.0);
    let serial =
        lddt_with_options(&model, &reference, &options, &execution(1)).expect("serial score");

    for workers in [2, 3, 4, 8, 16] {
        let parallel = lddt_with_options(&model, &reference, &options, &execution(workers))
            .expect("parallel score");
        assert!(
            (parallel - serial).abs() < f64::EPSILON,
            "worker count {workers} changed the score from {serial} to {parallel}"
        );
    }
}
