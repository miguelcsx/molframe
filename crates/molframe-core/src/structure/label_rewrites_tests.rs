use super::super::fixture::sample;
use super::{LabelRewrites, Structure};
use crate::index::ChainIndex;

#[test]
fn distinct_labels_report_nothing_to_rewrite() {
    assert!(sample().label_rewrites().is_empty());
    assert_eq!(sample().label_rewrites(), LabelRewrites::default());
}

#[test]
fn a_renamed_chain_is_reported_with_both_labels() {
    let mut data = sample().data().clone();
    let Ok(repeat) = data.dictionary.intern("A") else {
        panic!("intern");
    };
    assert!(
        data.topology
            .chains
            .relabel(ChainIndex::new(1), repeat)
            .is_ok()
    );
    let rewrites = Structure::new(data).label_rewrites();
    assert_eq!(rewrites.chains, [("A".to_owned(), "A2".to_owned())]);
    assert_eq!(rewrites.residues, 0);
    assert!(!rewrites.is_empty());
}
