use super::upgma;

#[test]
fn a_single_taxon_is_a_lone_leaf() {
    let tree = upgma(&["A"], &[vec![0.0]]);
    let Some(tree) = tree else {
        panic!("one taxon is valid");
    };
    assert_eq!(tree.to_newick(), "A;");
}

#[test]
fn two_close_taxa_join_before_a_distant_one() {
    // A and B are 2 apart; both are 4 from C. UPGMA joins A,B at height 1, then
    // that cluster with C at height 2.
    let labels = ["A", "B", "C"];
    let distances = vec![
        vec![0.0, 2.0, 4.0],
        vec![2.0, 0.0, 4.0],
        vec![4.0, 4.0, 0.0],
    ];
    let Some(tree) = upgma(&labels, &distances) else {
        panic!("valid matrix");
    };
    assert_eq!(tree.to_newick(), "((A:1,B:1):1,C:2);");
}

#[test]
fn a_non_square_matrix_is_rejected() {
    assert!(upgma(&["A", "B"], &[vec![0.0, 1.0]]).is_none());
}

#[test]
fn an_empty_input_has_no_tree() {
    assert!(upgma(&[], &[]).is_none());
}
