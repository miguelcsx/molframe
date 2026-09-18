use super::neighbor_joining;

#[test]
fn additive_distances_recover_the_generating_tree() {
    // Distances from a tree with edges A-x=1, B-x=2, x-y=1, y-C=3, y-D=4.
    let labels = ["A", "B", "C", "D"];
    let distances = vec![
        vec![0.0, 3.0, 5.0, 6.0],
        vec![3.0, 0.0, 6.0, 7.0],
        vec![5.0, 6.0, 0.0, 7.0],
        vec![6.0, 7.0, 7.0, 0.0],
    ];
    let Some(tree) = neighbor_joining(&labels, &distances) else {
        panic!("valid matrix");
    };
    // A and B join first with their true branch lengths of 1 and 2.
    assert_eq!(tree.to_newick(), "(((A:1,B:2):1,C:3):0,D:4);");
}

#[test]
fn two_taxa_join_directly() {
    let Some(tree) = neighbor_joining(&["A", "B"], &[vec![0.0, 4.0], vec![4.0, 0.0]]) else {
        panic!("valid");
    };
    assert_eq!(tree.to_newick(), "(A:0,B:4);");
}

#[test]
fn a_non_square_matrix_is_rejected() {
    assert!(neighbor_joining(&["A", "B"], &[vec![0.0]]).is_none());
}
