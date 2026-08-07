use super::*;

#[test]
fn a_list_reuses_candidates_only_within_half_the_skin() {
    let positions = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
    let list = match NeighborList::build(
        &positions,
        &[0, 1],
        &[0, 1],
        1.0,
        0.5,
        None,
        CoordinateGeneration::INITIAL,
    ) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    let moved = [[0.2, 0.0, 0.0], [1.2, 0.0, 0.0]];
    let pairs = match list.pairs(&moved, 1.0, None) {
        Ok(pairs) => pairs,
        Err(error) => panic!("reuse failed: {error}"),
    };
    assert_eq!(pairs.len(), 1);

    let stale = [[0.3, 0.0, 0.0], [1.2, 0.0, 0.0]];
    assert_eq!(
        list.pairs(&stale, 1.0, None),
        Err(SpatialError::StaleNeighborList)
    );
}

#[test]
fn generation_identity_is_available_to_structure_caches() {
    let positions = [[0.0, 0.0, 0.0]];
    let list = match NeighborList::build(
        &positions,
        &[0],
        &[0],
        1.0,
        0.5,
        None,
        CoordinateGeneration::INITIAL,
    ) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    assert!(list.is_current(CoordinateGeneration::INITIAL));
    assert!(!list.is_current(CoordinateGeneration::INITIAL.next()));
}
