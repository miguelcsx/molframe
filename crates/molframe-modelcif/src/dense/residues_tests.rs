use super::ResidueIndex;

#[test]
fn residues_keep_order_and_intern_labels() {
    let mut index = ResidueIndex::new();
    assert_eq!(index.push("A", 1, Some("MET")), 0);
    assert_eq!(index.push("A", 2, None), 1);
    assert_eq!(index.push("B", 1, Some("GLY")), 2);
    assert_eq!(index.len(), 3);
    let second = index.get(1).expect("residue");
    assert_eq!(
        (second.chain, second.sequence, second.component),
        ("A", 2, None)
    );
    assert_eq!(
        index.get(2).and_then(|residue| residue.component),
        Some("GLY")
    );
    assert_eq!(index.position("B", 1), Some(2));
    assert_eq!(index.position("B", 2), None);
    assert!(index.has_components());
    assert!(index.get(3).is_none());
}
