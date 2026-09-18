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

#[test]
fn cell_grid_coarsening_obeys_the_public_memory_policy() {
    let options = CellGridOptions {
        maximum_cell_count: 1,
        edge_growth_factor: 4.0,
    };
    let Ok((_, dimensions, count)) =
        geometry::grid_geometry([0.0, 0.0, 0.0], [3.0, 3.0, 3.0], 1.0, options)
    else {
        panic!("finite test geometry must be representable");
    };
    assert_eq!(dimensions, [1, 1, 1]);
    assert_eq!(count, 1);
}

#[test]
fn query_coordinates_below_target_bounds_clamp_to_the_first_cell() {
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let list = match CellList::build(&positions, &[2], 1.0, None) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    let pairs = match list.pairs(&[0, 1], 1.0) {
        Ok(pairs) => pairs,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(pairs.len(), 2);
}
