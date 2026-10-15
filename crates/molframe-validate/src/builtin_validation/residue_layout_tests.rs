use super::*;
use molframe_bench::{Sample, structure as sample};

#[test]
fn crambin_layout_has_forty_six_linked_residues() {
    let structure = sample(Sample::Tiny);
    let layout = Layout::new(&structure);
    assert_eq!(layout.len(), 46);
    assert!((0..46).all(|i| layout.is_amino_acid(i)));
    assert!(!layout.peptide_linked(0));
    assert!((1..46).all(|i| layout.peptide_linked(i)));
    assert!(!layout.is_nucleotide(3));
    assert!(layout.position(0, "CA").is_some());
    assert!(layout.position(0, "XX").is_none());
}

#[test]
fn chain_boundaries_break_the_link() {
    let structure = sample(Sample::Medium);
    let layout = Layout::new(&structure);
    let unlinked = (1..layout.len())
        .filter(|i| layout.is_amino_acid(*i) && !layout.peptide_linked(*i))
        .count();
    assert!(unlinked >= 3, "four chains give at least three starts");
}
