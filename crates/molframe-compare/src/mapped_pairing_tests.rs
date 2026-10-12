//! Naming the chain pairs to score, and natives with spare identical chains.

use super::tests::{
    HOMODIMER_NATIVE, HOMODIMER_SWAPPED, NATIVE, RENAMED, dockq_options, mapping, provider,
    structure,
};
use super::{mapped_dockq, mapped_dockq_for_chains};

#[test]
fn named_chain_pairs_are_scored_even_when_sequences_would_choose_otherwise() {
    let native = structure(HOMODIMER_NATIVE);
    let model = structure(HOMODIMER_SWAPPED);
    let provider = provider();
    let options = mapping(&provider);
    let (fits, comparison) = mapped_dockq_for_chains(
        &model,
        &native,
        ("A", "Y"),
        ("B", "X"),
        &options,
        dockq_options(),
    )
    .expect("the named pairs correspond");
    assert!((fits.score - 1.0).abs() < 1e-6, "score {}", fits.score);
    assert!(fits.fnonnat.abs() < 1e-12);
    assert_eq!(comparison.chain_assignments_tried, 1);
    let (forced, _) = mapped_dockq_for_chains(
        &model,
        &native,
        ("A", "X"),
        ("B", "Y"),
        &options,
        dockq_options(),
    )
    .expect("the named pairs correspond");
    assert!(forced.score < 0.9, "score {}", forced.score);
}

#[test]
fn a_named_pair_below_the_identity_floor_is_refused() {
    let native = structure(NATIVE);
    let model = structure(RENAMED);
    let provider = provider();
    let result = mapped_dockq_for_chains(
        &model,
        &native,
        ("A", "X"),
        ("B", "Y"),
        &mapping(&provider),
        dockq_options(),
    );
    assert!(matches!(result, Err(super::MappedCompareError::Mapping(_))));
}

#[test]
fn a_native_with_more_identical_chains_than_the_model_still_scores() {
    // The native holds three identical chains, the interface pair being the last
    // two; the model holds only the pair, so one native chain has no partner.
    let native = structure(
        "ATOM 1 C CA ALA A 1 50 50 50\n\
ATOM 2 C CA GLY A 2 52 50 50\n\
ATOM 3 C CA SER A 3 50 52 50\n\
ATOM 4 C CA ALA B 1 0 0 0\n\
ATOM 5 C CA GLY B 2 2 0 0\n\
ATOM 6 C CA SER B 3 0 2 0\n\
ATOM 7 C CA ALA C 1 1 1 3\n\
ATOM 8 C CA GLY C 2 3 1 3\n\
ATOM 9 C CA SER C 3 1 3 3\n",
    );
    let model = structure(HOMODIMER_SWAPPED);
    let provider = provider();
    let (score, _) = mapped_dockq(
        &model,
        &native,
        "B",
        "C",
        &mapping(&provider),
        dockq_options(),
    )
    .expect("the pair corresponds");
    assert!((score.score - 1.0).abs() < 1e-6, "score {}", score.score);
}
