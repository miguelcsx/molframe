use super::{condensed_pair, tm_score_all_pairs, tm_score_pairs};
use crate::failure::CompareError;
use crate::tm_score;
use molframe_core::ExecutionContext;

/// A deterministic wandering chain, perturbed per structure.
fn structure(seed: u32, atoms: usize) -> Vec<[f32; 3]> {
    let mut state = seed.wrapping_mul(2_654_435_761).wrapping_add(12_345);
    let mut next = move || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        f32::from(u16::try_from(state >> 16).unwrap_or(0)) / 65_536.0 - 0.5
    };
    let mut point = [0.0_f32; 3];
    (0..atoms)
        .map(|_| {
            point = [
                point[0] + 3.0 * next() + 1.5,
                point[1] + 3.0 * next(),
                point[2] + 3.0 * next(),
            ];
            point
        })
        .collect()
}

fn context(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("a valid worker budget")
}

#[test]
fn condensed_index_round_trips_for_every_pair() {
    for count in 2..40 {
        let mut expected = Vec::new();
        for row in 0..count {
            for column in row + 1..count {
                expected.push((row, column));
            }
        }
        for (index, pair) in expected.iter().enumerate() {
            assert_eq!(
                condensed_pair(count, index),
                Some(*pair),
                "n={count} k={index}"
            );
        }
        assert_eq!(condensed_pair(count, expected.len()), None);
    }
}

#[test]
fn batch_scores_equal_single_pair_scores_at_every_worker_count() {
    let owned: Vec<Vec<[f32; 3]>> = (0..6).map(|seed| structure(seed, 40)).collect();
    let views: Vec<&[[f32; 3]]> = owned.iter().map(Vec::as_slice).collect();
    let mut expected = Vec::new();
    for row in 0..views.len() {
        for column in row + 1..views.len() {
            expected.push(
                tm_score(views[row], views[column])
                    .expect("scorable pair")
                    .to_bits(),
            );
        }
    }
    for workers in [1, 2, 5, 16] {
        let found = tm_score_all_pairs(&views, &context(workers)).expect("batch succeeds");
        let bits: Vec<u64> = found.iter().map(|value| value.to_bits()).collect();
        assert_eq!(bits, expected, "{workers} workers");
    }
}

#[test]
fn explicit_pairs_follow_the_requested_order_and_direction() {
    let owned: Vec<Vec<[f32; 3]>> = (0..3).map(|seed| structure(seed, 30)).collect();
    let views: Vec<&[[f32; 3]]> = owned.iter().map(Vec::as_slice).collect();
    let pairs = [(2, 0), (0, 2), (1, 1)];
    let found = tm_score_pairs(&views, &pairs, &context(4)).expect("batch succeeds");
    for (value, &(a, b)) in found.iter().zip(&pairs) {
        let single = tm_score(views[a], views[b]).expect("scorable pair");
        assert_eq!(value.to_bits(), single.to_bits());
    }
}

#[test]
fn bad_indices_and_unequal_lengths_are_reported() {
    let short = structure(1, 20);
    let long = structure(2, 30);
    let views: Vec<&[[f32; 3]]> = vec![&short, &long];
    assert!(matches!(
        tm_score_pairs(&views, &[(0, 2)], &context(2)),
        Err(CompareError::InvalidMapping)
    ));
    assert!(matches!(
        tm_score_all_pairs(&views, &context(2)),
        Err(CompareError::LengthMismatch { .. })
    ));
    assert!(
        tm_score_all_pairs(&[], &context(2))
            .expect("empty input is valid")
            .is_empty()
    );
}

#[test]
fn sampled_searches_above_the_seed_cap_also_match_single_pair_scores() {
    let owned: Vec<Vec<[f32; 3]>> = (0..3).map(|seed| structure(seed, 1_100)).collect();
    let views: Vec<&[[f32; 3]]> = owned.iter().map(Vec::as_slice).collect();
    let found = tm_score_all_pairs(&views, &context(4)).expect("batch succeeds");
    let single = [(0, 1), (0, 2), (1, 2)].map(|(a, b)| {
        tm_score(views[a], views[b])
            .expect("scorable pair")
            .to_bits()
    });
    let bits: Vec<u64> = found.iter().map(|value| value.to_bits()).collect();
    assert_eq!(bits, single);
}
