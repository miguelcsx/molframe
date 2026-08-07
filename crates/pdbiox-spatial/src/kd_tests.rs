use super::*;

#[test]
fn nearest_neighbours_are_sorted_by_distance_then_atom() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ];
    let tree = match KdTree::build(&positions, &[0, 1, 2, 3], None) {
        Ok(tree) => tree,
        Err(error) => panic!("build failed: {error}"),
    };
    let nearest = match tree.k_nearest(0, 2) {
        Ok(nearest) => nearest,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(
        nearest.iter().map(|item| item.0).collect::<Vec<_>>(),
        vec![1, 2]
    );
}

#[test]
fn radius_search_excludes_self_and_deduplicates_overlapping_groups() {
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
    let tree = match KdTree::build(&positions, &[0, 1], None) {
        Ok(tree) => tree,
        Err(error) => panic!("build failed: {error}"),
    };
    let pairs = match tree.pairs(&[0, 1], 1.0) {
        Ok(pairs) => pairs,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(pairs.len(), 1);
    assert_eq!((pairs[0].first, pairs[0].second), (0, 1));
}
