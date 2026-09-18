use super::{LadderDirection, NewickError, RerootError, TraversalOrder, Tree};

#[test]
fn a_leaf_renders_as_its_name() {
    let tree = Tree::Leaf {
        name: "A".to_string(),
    };
    assert_eq!(tree.to_newick(), "A;");
}

#[test]
fn a_clade_renders_with_branch_lengths() {
    let tree = Tree::Clade {
        left: Box::new(Tree::Leaf {
            name: "A".to_string(),
        }),
        left_length: 1.0,
        right: Box::new(Tree::Leaf {
            name: "B".to_string(),
        }),
        right_length: 1.5,
    };
    assert_eq!(tree.to_newick(), "(A:1,B:1.5);");
}

#[test]
fn a_written_tree_parses_back_to_itself() {
    let newick = "((A:1,B:1):1,C:2);";
    let tree = match Tree::from_newick(newick) {
        Ok(tree) => tree,
        Err(error) => panic!("valid Newick: {error:?}"),
    };
    assert_eq!(tree.to_newick(), newick);
}

#[test]
fn a_bare_leaf_parses() {
    let tree = match Tree::from_newick("A;") {
        Ok(tree) => tree,
        Err(error) => panic!("valid: {error:?}"),
    };
    assert_eq!(
        tree,
        Tree::Leaf {
            name: "A".to_string()
        }
    );
}

#[test]
fn a_truncated_tree_is_an_error() {
    assert!(matches!(
        Tree::from_newick("((A:1,B:1)"),
        Err(NewickError::UnexpectedEnd | NewickError::Unexpected(_))
    ));
}

#[test]
fn traversal_orders_visit_every_node_in_the_requested_order() {
    let Ok(tree) = Tree::from_newick("((A:1,B:1):2,C:3);") else {
        panic!("valid tree");
    };
    let preorder = tree.traverse(TraversalOrder::Preorder);
    let postorder = tree.traverse(TraversalOrder::Postorder);
    assert_eq!(preorder.len(), 5);
    assert!(matches!(preorder[0], Tree::Clade { .. }));
    assert!(matches!(postorder[4], Tree::Clade { .. }));
}

#[test]
fn ladderisation_is_deterministic_for_equal_and_unequal_clades() {
    let Ok(mut tree) = Tree::from_newick("(C:1,(B:1,A:1):1);") else {
        panic!("valid tree");
    };
    tree.ladderize(LadderDirection::Descending);
    assert_eq!(tree.to_newick(), "((A:1,B:1):1,C:1);");
}

#[test]
fn rerooting_preserves_leaf_distances_and_splits_the_selected_edge() {
    let Ok(tree) = Tree::from_newick("((A:1,B:2):3,C:4);") else {
        panic!("valid tree");
    };
    let Ok(rerooted) = tree.reroot_at_leaf("B", 0.25) else {
        panic!("valid reroot");
    };
    let Tree::Clade {
        left_length,
        right_length,
        ..
    } = rerooted
    else {
        panic!("reroot creates a clade");
    };
    assert!((left_length - 0.5).abs() < f64::EPSILON);
    assert!((right_length - 1.5).abs() < f64::EPSILON);
}

#[test]
fn rerooting_rejects_missing_and_ambiguous_leaf_names() {
    let Ok(missing) = Tree::from_newick("(A:1,B:1);") else {
        panic!("valid tree");
    };
    assert_eq!(
        missing.reroot_at_leaf("C", 0.5),
        Err(RerootError::LeafNotFound)
    );
    let Ok(duplicate) = Tree::from_newick("(A:1,A:1);") else {
        panic!("valid tree");
    };
    assert_eq!(
        duplicate.reroot_at_leaf("A", 0.5),
        Err(RerootError::AmbiguousLeaf)
    );
}
