use super::*;
use crate::Scoring;
use crate::align::dynamic::{Mode, run};

#[test]
fn linear_space_global_matches_full_affine_scores_over_generated_inputs() {
    let scoring = Scoring {
        match_score: 2,
        mismatch_score: -1,
        gap_open: -3,
        gap_extend: -1,
    };
    let alphabet = [b'A', b'B'];
    for left_len in 0..=6 {
        for right_len in 0..=6 {
            for left_bits in 0..(1usize << left_len) {
                for right_bits in 0..(1usize << right_len) {
                    let left = sequence(left_bits, left_len, alphabet);
                    let right = sequence(right_bits, right_len, alphabet);
                    let Ok(linear) = global_linear(
                        &left,
                        &right,
                        &scoring,
                        scoring.gap_open,
                        scoring.gap_extend,
                    ) else {
                        panic!("generated linear alignment overflowed")
                    };
                    let Ok(full) = run(
                        &left,
                        &right,
                        &scoring,
                        scoring.gap_open,
                        scoring.gap_extend,
                        Mode::Global,
                        None,
                    ) else {
                        panic!("generated dynamic alignment overflowed")
                    };
                    assert_eq!(linear.score, full.score, "left={left:?} right={right:?}");
                    assert_eq!(
                        score_columns(&left, &right, &linear.columns, scoring),
                        linear.score
                    );
                }
            }
        }
    }
}

#[test]
fn workspace_is_bounded_by_the_shorter_axis() {
    let scoring = Scoring::simple();
    let long = vec![b'A'; 20_000];
    let short = vec![b'A'; 7];
    let Ok(alignment) = global_linear(
        &long,
        &short,
        &scoring,
        scoring.gap_open,
        scoring.gap_extend,
    ) else {
        panic!("bounded linear alignment overflowed")
    };
    assert_eq!(alignment.columns.len(), long.len());
}

fn sequence(bits: usize, length: usize, alphabet: [u8; 2]) -> Vec<u8> {
    (0..length)
        .map(|index| alphabet[(bits >> index) & 1])
        .collect()
}

fn score_columns(left: &[u8], right: &[u8], columns: &[Column], scoring: Scoring) -> i32 {
    let mut score = 0;
    let mut gap = None;
    for column in columns {
        match (column.left, column.right) {
            (Some(first), Some(second)) => {
                score += scoring.score(left[first], right[second]);
                gap = None;
            }
            (Some(_), None) => {
                score += if gap == Some(State::RightGap) {
                    scoring.gap_extend
                } else {
                    scoring.gap_open
                };
                gap = Some(State::RightGap);
            }
            (None, Some(_)) => {
                score += if gap == Some(State::LeftGap) {
                    scoring.gap_extend
                } else {
                    scoring.gap_open
                };
                gap = Some(State::LeftGap);
            }
            (None, None) => panic!("alignment column cannot contain two gaps"),
        }
    }
    score
}
