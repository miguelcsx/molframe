use super::Scoring;

#[test]
fn equal_residues_score_the_match_value() {
    let scoring = Scoring::simple();
    assert_eq!(scoring.substitution(b'A', b'A'), 1);
}

#[test]
fn different_residues_score_the_mismatch_value() {
    let scoring = Scoring::simple();
    assert_eq!(scoring.substitution(b'A', b'C'), -1);
}
