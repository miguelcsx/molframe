use super::*;
use crate::Scoring;

#[test]
fn similar_words_are_score_then_lexicographically_ranked() {
    let options = SimilarKmerOptions {
        minimum_score: 0,
        max_results: 4,
        max_candidates: 16,
    };
    let Ok(words) = similar_kmers(b"AA", b"AC", &Scoring::simple(), options) else {
        panic!("bounded enumeration");
    };
    assert_eq!(
        words[0],
        SimilarKmer {
            word: b"AA".to_vec(),
            score: 2
        }
    );
    assert_eq!(words.len(), 3);
}

#[test]
fn an_incomplete_cartesian_search_is_refused() {
    let options = SimilarKmerOptions {
        minimum_score: 0,
        max_results: 1,
        max_candidates: 3,
    };
    assert_eq!(
        similar_kmers(b"AA", b"AC", &Scoring::simple(), options),
        Err(SimilarKmerError::CandidateLimit {
            required: 4,
            allowed: 3
        })
    );
}

#[test]
fn an_unrepresentable_candidate_score_is_refused() {
    let scoring = Scoring {
        match_score: i32::MAX,
        mismatch_score: 0,
        gap_open: 0,
        gap_extend: 0,
    };
    let options = SimilarKmerOptions {
        minimum_score: 0,
        max_results: 1,
        max_candidates: 1,
    };
    assert_eq!(
        similar_kmers(b"AA", b"A", &scoring, options),
        Err(SimilarKmerError::NumericOverflow)
    );
}
