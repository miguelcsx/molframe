use super::*;

#[test]
fn bit_parallel_lcs_matches_scalar_reference_across_word_boundaries() {
    for left_length in [0, 1, 7, 63, 64, 65, 127, 129] {
        for right_length in [0, 2, 31, 64, 66, 128] {
            let left = sequence(left_length, 3);
            let mut right = sequence(right_length, 5);
            for position in (11..right.len()).step_by(17) {
                right[position] = b'Y';
            }
            let maximum = left.len().max(right.len());
            let sequences = [left.as_slice(), right.as_slice()];
            let (slots, alphabet_size) = symbol_slots(&sequences);
            let mut workspace =
                LcsWorkspace::new(maximum, &slots, alphabet_size, usize::MAX, usize::MAX)
                    .unwrap_or_else(|error| panic!("workspace failed: {error}"));
            let actual = workspace
                .distance(&left, &right)
                .unwrap_or_else(|error| panic!("LCS failed: {error}"));
            let expected = scalar_distance(&left, &right);
            assert!(
                (actual - expected).abs() < f64::EPSILON,
                "{left_length}x{right_length}"
            );
        }
    }
}

#[test]
fn gaps_do_not_contribute_to_guide_distance() {
    let sequences = [b"A-CGT".as_slice(), b"ACG-T".as_slice()];
    let (slots, alphabet_size) = symbol_slots(&sequences);
    let mut workspace = LcsWorkspace::new(8, &slots, alphabet_size, usize::MAX, usize::MAX)
        .unwrap_or_else(|error| panic!("workspace failed: {error}"));
    assert_eq!(workspace.distance(b"A-CGT", b"ACG-T"), Ok(0.0));
}

#[test]
fn lcs_masks_scale_with_the_observed_alphabet() {
    let sequences = [b"ACDEFGHIKLMNPQRSTVWY".as_slice()];
    let (slots, alphabet_size) = symbol_slots(&sequences);
    let workspace = LcsWorkspace::new(1_000_000, &slots, alphabet_size, usize::MAX, usize::MAX)
        .unwrap_or_else(|error| panic!("workspace failed: {error}"));
    let words = 1_000_000_usize.div_ceil(64);
    assert_eq!(alphabet_size, 20);
    assert_eq!(workspace.masks.len(), words * alphabet_size);
    assert!(workspace.masks.len() * 10 < words * 256);
}

#[test]
fn mst_edges_produce_a_deterministic_single_linkage_tree() {
    let edges = [
        MstEdge {
            distance: 0.1,
            first: 0,
            second: 1,
        },
        MstEdge {
            distance: 0.2,
            first: 1,
            second: 2,
        },
        MstEdge {
            distance: 0.3,
            first: 2,
            second: 3,
        },
    ];
    assert_eq!(merges_from_mst(&edges, 4), [(0, 1), (0, 2), (0, 3)]);
}

fn sequence(length: usize, stride: usize) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    (0..length)
        .map(|index| ALPHABET[(index * stride + index / 7) % ALPHABET.len()])
        .collect()
}

fn scalar_distance(left: &[u8], right: &[u8]) -> f64 {
    let left = left
        .iter()
        .copied()
        .filter(|symbol| *symbol != GAP)
        .collect::<Vec<_>>();
    let right = right
        .iter()
        .copied()
        .filter(|symbol| *symbol != GAP)
        .collect::<Vec<_>>();
    let mut previous = vec![0_u32; right.len() + 1];
    let mut current = previous.clone();
    for first in &left {
        for (column, second) in right.iter().enumerate() {
            current[column + 1] = if first == second {
                previous[column] + 1
            } else {
                previous[column + 1].max(current[column])
            };
        }
        std::mem::swap(&mut previous, &mut current);
        current.fill(0);
    }
    let normalizer = left.len().max(right.len());
    if normalizer == 0 {
        0.0
    } else {
        let normalizer = u32::try_from(normalizer)
            .unwrap_or_else(|error| panic!("test normalizer failed: {error}"));
        1.0 - f64::from(previous[right.len()]) / f64::from(normalizer)
    }
}
