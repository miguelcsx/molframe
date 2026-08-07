use super::*;

#[test]
fn contiguous_cells_find_only_pairs_inside_the_cutoff() {
    let positions = [[0.0, 0.0, 0.0], [0.9, 0.0, 0.0], [4.0, 0.0, 0.0]];
    let list = match CellList::build(&positions, &[0, 1, 2], 1.0, None) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    let pairs = match list.pairs(&[0, 1, 2], 1.0) {
        Ok(pairs) => pairs,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(pairs.len(), 1);
    assert_eq!((pairs[0].first, pairs[0].second), (0, 1));
}

#[test]
fn a_query_cannot_exceed_the_cutoff_the_grid_was_built_for() {
    let positions = [[0.0, 0.0, 0.0]];
    let list = match CellList::build(&positions, &[0], 1.0, None) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    assert_eq!(list.pairs(&[0], 2.0), Err(SpatialError::InvalidCutoff));
}
