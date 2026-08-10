use super::{ORDER, blosum62};
use crate::scoring::Score;

#[test]
fn a_known_diagonal_entry_matches_the_published_matrix() {
    let m = blosum62();
    assert_eq!(m.get(b'W', b'W'), 11);
    assert_eq!(m.get(b'A', b'A'), 4);
    assert_eq!(m.get(b'C', b'C'), 9);
}

#[test]
fn a_known_off_diagonal_entry_matches_the_published_matrix() {
    let m = blosum62();
    assert_eq!(m.get(b'A', b'W'), -3);
    assert_eq!(m.get(b'F', b'Y'), 3);
    assert_eq!(m.get(b'D', b'E'), 2);
}

#[test]
fn the_matrix_is_symmetric_across_the_alphabet() {
    let m = blosum62();
    for &a in &ORDER {
        for &b in &ORDER {
            assert_eq!(m.get(a, b), m.get(b, a));
        }
    }
}

#[test]
fn lower_case_input_scores_the_same_as_upper_case() {
    let m = blosum62();
    assert_eq!(m.get(b'w', b'w'), 11);
    assert_eq!(m.get(b'a', b'W'), -3);
}

#[test]
fn a_byte_outside_the_alphabet_is_scored_as_any() {
    let m = blosum62();
    // '?' is not a residue, so it falls back to the 'X' column.
    assert_eq!(m.get(b'?', b'A'), m.get(b'X', b'A'));
}

#[test]
fn the_score_trait_agrees_with_the_inherent_lookup() {
    let m = blosum62();
    assert_eq!(Score::score(&m, b'H', b'Y'), m.get(b'H', b'Y'));
}
