use super::{EmptyLddtPolicy, LddtOptions, lddt, lddt_with_options};
use crate::failure::CompareError;

#[test]
fn an_identical_model_reproduces_every_distance() {
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
    let score = match lddt(&reference, &reference, 15.0) {
        Ok(score) => score,
        Err(error) => panic!("valid: {error}"),
    };
    assert!((score - 1.0).abs() < 1e-12);
}

#[test]
fn a_distance_broken_past_every_tolerance_scores_zero() {
    let reference = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [13.0, 0.0, 0.0]];
    let score = match lddt(&model, &reference, 15.0) {
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
    let score = match lddt(&model, &reference, 15.0) {
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
    let score = match lddt(&model, &reference, 2.0) {
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
        lddt(&model, &reference, 15.0),
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
            }
        ),
        Err(CompareError::NoComparablePairs)
    ));
}
