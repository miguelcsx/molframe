use super::seed_and_extend;
use crate::scoring::Scoring;

#[test]
fn a_shared_word_seeds_a_local_alignment() {
    let Ok(Some(alignment)) = seed_and_extend(b"CGT", b"AACGTT", 2, Scoring::simple()) else {
        panic!("CG and GT are shared seeds");
    };
    assert_eq!(alignment.score, 3);
}

#[test]
fn no_shared_word_means_no_hit() {
    assert!(seed_and_extend(b"TTTT", b"AAAA", 2, Scoring::simple()).is_ok_and(|hit| hit.is_none()));
}

#[test]
fn a_seed_longer_than_the_query_cannot_match() {
    assert!(seed_and_extend(b"AC", b"AACGTT", 3, Scoring::simple()).is_ok_and(|hit| hit.is_none()));
}
