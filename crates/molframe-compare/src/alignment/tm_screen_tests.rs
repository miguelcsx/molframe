use super::{TmPrefilter, tm_screen_all_pairs, tm_screen_pairs};
use crate::failure::CompareError;
use crate::tm_bound::tests::{chain, nudged};
use crate::{tm_score_all_pairs, tm_score_pairs};
use molframe_core::ExecutionContext;

fn context(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("a valid worker budget")
}

/// Three families of near-copies plus unrelated singletons.
fn library() -> Vec<Vec<[f32; 3]>> {
    let mut out = Vec::new();
    for family in 0..3_u32 {
        let base = chain(50 + family, 60);
        for member in 0..3_u32 {
            out.push(nudged(&base, family * 10 + member, 0.4));
        }
    }
    for seed in 0..4_u32 {
        out.push(chain(200 + seed, 60));
    }
    out
}

fn views(owned: &[Vec<[f32; 3]>]) -> Vec<&[[f32; 3]]> {
    owned.iter().map(Vec::as_slice).collect()
}

#[test]
fn exact_screen_keeps_every_pair_at_or_above_the_threshold() {
    let owned = library();
    let views = views(&owned);
    let full = tm_score_all_pairs(&views, &context(2)).expect("full scores");
    for threshold in [0.3, 0.5, 0.7, 0.9] {
        let filter = TmPrefilter::new(threshold).expect("valid");
        let screened = tm_screen_all_pairs(&views, &filter, &context(2)).expect("screened");
        assert_eq!(screened.scores.len(), full.len());
        let mut kept = 0;
        for (index, (&exact, &got)) in full.iter().zip(&screened.scores).enumerate() {
            if exact >= threshold {
                assert_eq!(
                    got.to_bits(),
                    exact.to_bits(),
                    "pair {index} at {threshold}"
                );
                kept += 1;
            } else if !got.is_nan() {
                assert_eq!(got.to_bits(), exact.to_bits());
            }
        }
        assert!(kept > 0, "fixture must have hits at {threshold}");
        let total = u64::try_from(full.len()).expect("fits");
        assert_eq!(
            screened.evaluated + screened.pruned_bound + screened.pruned_coarse,
            total
        );
        assert_eq!(screened.pruned_coarse, 0);
    }
}

#[test]
fn exact_screen_prunes_unrelated_pairs() {
    let owned = library();
    let views = views(&owned);
    let filter = TmPrefilter::new(0.9).expect("valid");
    let screened = tm_screen_all_pairs(&views, &filter, &context(1)).expect("screened");
    assert!(
        screened.pruned_bound > 0,
        "the bound should prune something"
    );
    assert!(screened.scores.iter().any(|v| v.is_nan()));
}

#[test]
fn screen_is_identical_at_every_worker_count() {
    let owned = library();
    let views = views(&owned);
    let filter = TmPrefilter::new(0.5)
        .and_then(|f| f.with_coarse_ratio(0.5))
        .expect("valid");
    let baseline = tm_screen_all_pairs(&views, &filter, &context(1)).expect("screened");
    for workers in [2, 5] {
        let other = tm_screen_all_pairs(&views, &filter, &context(workers)).expect("screened");
        let same = baseline
            .scores
            .iter()
            .zip(&other.scores)
            .all(|(a, b)| a.to_bits() == b.to_bits());
        assert!(same, "workers={workers}");
        assert_eq!(baseline.evaluated, other.evaluated);
        assert_eq!(baseline.pruned_bound, other.pruned_bound);
        assert_eq!(baseline.pruned_coarse, other.pruned_coarse);
    }
}

#[test]
fn coarse_stage_only_prunes_pairs_the_single_fit_scores_low() {
    let owned = library();
    let views = views(&owned);
    let exact = TmPrefilter::new(0.5).expect("valid");
    let coarse = exact.with_coarse_ratio(1.0).expect("valid");
    let before = tm_screen_all_pairs(&views, &exact, &context(1)).expect("screened");
    let after = tm_screen_all_pairs(&views, &coarse, &context(1)).expect("screened");
    assert!(after.evaluated <= before.evaluated);
    assert_eq!(
        after.evaluated + after.pruned_bound + after.pruned_coarse,
        before.evaluated + before.pruned_bound
    );
    for (b, a) in before.scores.iter().zip(&after.scores) {
        if !a.is_nan() {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}

#[test]
fn pair_list_matches_the_condensed_screen() {
    let owned = library();
    let views = views(&owned);
    let filter = TmPrefilter::new(0.6).expect("valid");
    let pairs = [(0, 1), (0, 4), (3, 9), (5, 11)];
    let screened = tm_screen_pairs(&views, &pairs, &filter, &context(2)).expect("screened");
    let full = tm_score_pairs(&views, &pairs, &context(2)).expect("full");
    for (exact, got) in full.iter().zip(&screened.scores) {
        if *exact >= 0.6 {
            assert_eq!(got.to_bits(), exact.to_bits());
        }
    }
}

#[test]
fn invalid_inputs_are_rejected() {
    for bad in [0.0, -0.1, 1.5, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            TmPrefilter::new(bad),
            Err(CompareError::InvalidScoreInput)
        ));
    }
    let filter = TmPrefilter::new(0.5).expect("valid");
    assert!(matches!(
        filter.with_coarse_ratio(0.0),
        Err(CompareError::InvalidScoreInput)
    ));
    let owned = library();
    let views = views(&owned);
    assert!(matches!(
        tm_screen_pairs(&views, &[(0, 99)], &filter, &context(1)),
        Err(CompareError::InvalidMapping)
    ));
    let short = chain(1, 20);
    let long = chain(1, 30);
    let mixed: Vec<&[[f32; 3]]> = vec![&short, &long];
    assert!(matches!(
        tm_screen_all_pairs(&mixed, &filter, &context(1)),
        Err(CompareError::LengthMismatch { .. })
    ));
}

#[test]
fn coarse_stage_keeps_every_hit_on_family_libraries() {
    // Recall of the coarse stage is a property of the data; this library has
    // global matches, where one fit is enough.
    let owned = library();
    let views = views(&owned);
    let full = tm_score_all_pairs(&views, &context(2)).expect("full scores");
    let filter = TmPrefilter::new(0.5)
        .and_then(|f| f.with_coarse_ratio(0.5))
        .expect("valid");
    let screened = tm_screen_all_pairs(&views, &filter, &context(2)).expect("screened");
    for (exact, got) in full.iter().zip(&screened.scores) {
        if *exact >= 0.5 {
            assert_eq!(got.to_bits(), exact.to_bits());
        }
    }
}
