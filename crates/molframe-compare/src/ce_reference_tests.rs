//! Matrix-backed CE oracle retained only for bitwise differential tests.

use super::search::{Candidate, Fragment, find_paths};
use super::{CeOptions, CeSignificanceProfile, workspace};
use crate::numeric::usize_to_f64;

#[test]
fn banded_search_is_bitwise_equivalent_to_the_matrix_algorithm() {
    for window_size in 3..=8 {
        let seed = u16::try_from(window_size)
            .unwrap_or_else(|error| panic!("window does not fit the fixture seed: {error}"));
        let reference = coordinates(36 + window_size, seed);
        let mut mobile: Vec<_> = reference
            .iter()
            .enumerate()
            .map(|(index, point)| {
                let phase = u8::try_from(index % 5)
                    .unwrap_or_else(|error| panic!("fixture phase is invalid: {error}"));
                let drift = f32::from(phase) * 0.003;
                [-point[1] + 4.0 + drift, point[0] - 2.0, point[2] + 1.0]
            })
            .collect();
        mobile.insert(window_size + 2, [73.0, -41.0, 19.0]);
        let mut options = CeOptions {
            window_size,
            max_gap: 4,
            max_paths: 7,
            fragment_similarity_threshold: -3.0,
            path_similarity_threshold: -4.0,
            significance: if window_size == 8 {
                Some(CeSignificanceProfile::OriginalWindowEight)
            } else {
                None
            },
            memory_limit_bytes: 500_000_000,
        };
        if window_size % 2 == 0 {
            let plan = workspace::memory_plan(reference.len(), mobile.len(), options)
                .unwrap_or_else(|error| panic!("CE memory planning failed: {error}"));
            options.memory_limit_bytes = plan.minimum_bytes + 4_096;
            let bounded = workspace::memory_plan(reference.len(), mobile.len(), options)
                .unwrap_or_else(|error| panic!("bounded CE memory planning failed: {error}"));
            assert_eq!(bounded.cache_bytes, None);
        }
        let actual = find_paths(&reference, &mobile, options)
            .unwrap_or_else(|error| panic!("banded CE search failed: {error}"));
        let expected = matrix_find_paths(&reference, &mobile, options);
        assert_eq!(actual.len(), expected.len(), "window={window_size}");
        for (actual, expected) in actual.iter().zip(&expected) {
            assert_eq!(actual.fragments, expected.fragments, "window={window_size}");
            assert_eq!(
                actual.similarity.to_bits(),
                expected.similarity.to_bits(),
                "window={window_size}"
            );
        }
    }
}

fn coordinates(count: usize, seed: u16) -> Vec<[f32; 3]> {
    (0..count)
        .map(|index| {
            let index = u16::try_from(index)
                .unwrap_or_else(|error| panic!("fixture index is too large: {error}"));
            let value = f32::from(index) + f32::from(seed) * 0.17;
            [
                value * 1.31,
                (value * 0.37).sin() * 3.0 + value * 0.02,
                (value * 0.23).cos() * 2.0 - value * 0.01,
            ]
        })
        .collect()
}

struct MatrixDistances {
    size: usize,
    values: Vec<f64>,
}

impl MatrixDistances {
    fn new(points: &[[f32; 3]]) -> Self {
        let mut values = vec![0.0; points.len() * points.len()];
        for row in 0..points.len() {
            for column in row + 1..points.len() {
                let distance = molframe_geom::distance(points[row], points[column]);
                values[row * points.len() + column] = distance;
                values[column * points.len() + row] = distance;
            }
        }
        Self {
            size: points.len(),
            values,
        }
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.size + column]
    }
}

struct MatrixSimilarities {
    columns: usize,
    values: Vec<f64>,
}

impl MatrixSimilarities {
    fn new(reference: &MatrixDistances, mobile: &MatrixDistances, window: usize) -> Self {
        let rows = reference.size - window + 1;
        let columns = mobile.size - window + 1;
        let mut values = Vec::with_capacity(rows * columns);
        let terms = (window - 1) * (window - 2) / 2;
        for reference_start in 0..rows {
            for mobile_start in 0..columns {
                let mut difference = 0.0;
                for left in 0..window - 2 {
                    for right in left + 2..window {
                        difference += (reference
                            .get(reference_start + left, reference_start + right)
                            - mobile.get(mobile_start + left, mobile_start + right))
                        .abs();
                    }
                }
                values.push(-difference / usize_to_f64(terms));
            }
        }
        Self { columns, values }
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.columns + column]
    }
}

fn matrix_find_paths(
    reference_points: &[[f32; 3]],
    mobile_points: &[[f32; 3]],
    options: CeOptions,
) -> Vec<Candidate> {
    let reference = MatrixDistances::new(reference_points);
    let mobile = MatrixDistances::new(mobile_points);
    let similarity = MatrixSimilarities::new(&reference, &mobile, options.window_size);
    let reference_starts = reference.size - options.window_size + 1;
    let mobile_starts = mobile.size - options.window_size + 1;
    let mut best = Vec::new();
    for reference_start in 0..reference_starts {
        if matrix_cannot_beat(
            reference_start,
            reference.size,
            &best,
            options.window_size,
            options.max_paths,
        ) {
            break;
        }
        for mobile_start in 0..mobile_starts {
            if matrix_cannot_beat(
                mobile_start,
                mobile.size,
                &best,
                options.window_size,
                options.max_paths,
            ) {
                break;
            }
            let initial = similarity.get(reference_start, mobile_start);
            if initial <= options.fragment_similarity_threshold {
                continue;
            }
            let mut candidate = Candidate {
                fragments: vec![Fragment {
                    reference: reference_start,
                    mobile: mobile_start,
                }],
                similarity: initial,
            };
            matrix_extend_path(&mut candidate, &similarity, &reference, &mobile, options);
            matrix_insert_candidate(&mut best, candidate, options.max_paths);
        }
    }
    best
}

fn matrix_cannot_beat(
    start: usize,
    length: usize,
    best: &[Candidate],
    window: usize,
    max_paths: usize,
) -> bool {
    if best.len() < max_paths {
        return false;
    }
    best.last().is_some_and(|worst| {
        let gaps = worst.fragments.len().saturating_sub(1);
        let Some(span) = window.checked_mul(gaps) else {
            return true;
        };
        let Some(last_start) = length.checked_sub(span) else {
            return true;
        };
        start > last_start
    })
}

#[test]
fn length_pruning_waits_until_the_top_k_is_full() {
    let candidate = Candidate {
        fragments: vec![
            Fragment {
                reference: 0,
                mobile: 0,
            };
            4
        ],
        similarity: -1.0,
    };
    assert!(!super::search::cannot_beat(9, 32, &[candidate], 8, 2));

    let full = Candidate {
        fragments: vec![
            Fragment {
                reference: 0,
                mobile: 0,
            };
            4
        ],
        similarity: -1.0,
    };
    assert!(super::search::cannot_beat(9, 32, &[full], 8, 1));
}

fn matrix_extend_path(
    candidate: &mut Candidate,
    similarity: &MatrixSimilarities,
    reference: &MatrixDistances,
    mobile: &MatrixDistances,
    options: CeOptions,
) {
    loop {
        let Some(last) = candidate.fragments.last().copied() else {
            return;
        };
        let mut best_extension: Option<(Fragment, f64)> = None;
        for gap_index in 0..=options.max_gap * 2 {
            let mut next = Fragment {
                reference: last.reference + options.window_size,
                mobile: last.mobile + options.window_size,
            };
            if (gap_index + 1).is_multiple_of(2) {
                next.reference += gap_index.div_ceil(2);
            } else {
                next.mobile += gap_index.div_ceil(2);
            }
            if next.reference + options.window_size > reference.size
                || next.mobile + options.window_size > mobile.size
                || similarity.get(next.reference, next.mobile)
                    <= options.fragment_similarity_threshold
            {
                continue;
            }
            let cross = candidate
                .fragments
                .iter()
                .map(|fragment| {
                    matrix_cross_similarity(*fragment, next, reference, mobile, options.window_size)
                })
                .sum::<f64>()
                / usize_to_f64(candidate.fragments.len());
            if cross > options.path_similarity_threshold
                && best_extension.is_none_or(|(_, best_cross)| cross > best_cross)
            {
                best_extension = Some((next, cross));
            }
        }
        let Some((next, cross)) = best_extension else {
            break;
        };
        let count = usize_to_f64(candidate.fragments.len());
        let current_terms = count + count * (count - 1.0) / 2.0;
        let new_terms = count + 1.0 + count * (count + 1.0) / 2.0;
        let next_similarity = (current_terms * candidate.similarity
            + count * cross
            + similarity.get(next.reference, next.mobile))
            / new_terms;
        if next_similarity <= options.path_similarity_threshold {
            break;
        }
        candidate.fragments.push(next);
        candidate.similarity = next_similarity;
    }
}

fn matrix_cross_similarity(
    first: Fragment,
    second: Fragment,
    reference: &MatrixDistances,
    mobile: &MatrixDistances,
    window: usize,
) -> f64 {
    let mut difference = (reference.get(first.reference, second.reference)
        - mobile.get(first.mobile, second.mobile))
    .abs();
    difference += (reference.get(first.reference + window - 1, second.reference + window - 1)
        - mobile.get(first.mobile + window - 1, second.mobile + window - 1))
    .abs();
    for offset in 1..window - 1 {
        difference += (reference.get(
            first.reference + offset,
            second.reference + window - 1 - offset,
        ) - mobile.get(first.mobile + offset, second.mobile + window - 1 - offset))
        .abs();
    }
    -difference / usize_to_f64(window)
}

fn matrix_insert_candidate(best: &mut Vec<Candidate>, candidate: Candidate, limit: usize) {
    best.push(candidate);
    best.sort_by(|left, right| {
        right
            .fragments
            .len()
            .cmp(&left.fragments.len())
            .then_with(|| right.similarity.total_cmp(&left.similarity))
    });
    best.truncate(limit);
}
