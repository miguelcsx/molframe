use super::*;

#[test]
fn cross_indices_remain_local_even_when_the_right_index_is_smaller() {
    let first = [[10.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let second = [[0.5, 0.0, 0.0]];
    let pairs = cross_pairs(
        &first,
        &second,
        1.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    )
    .expect("valid independent coordinate sets");
    assert_eq!(
        pairs,
        vec![CrossPair {
            first: 1,
            second: 0,
            distance_squared: 0.25
        }]
    );
}

#[test]
fn dense_sets_do_not_emit_within_set_pairs() {
    let first = [[0.0; 3]; 4];
    let second = [[0.0; 3]; 3];
    let pairs = cross_pairs(
        &first,
        &second,
        1.0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    )
    .expect("valid dense coordinate sets");
    assert_eq!(pairs.len(), first.len() * second.len());
    assert!(pairs.iter().all(|pair| pair.first < 4 && pair.second < 3));
}
