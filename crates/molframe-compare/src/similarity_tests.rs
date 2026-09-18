use super::contact_map_similarity;

#[test]
fn identical_maps_are_perfectly_similar() {
    let map = [(0u32, 1u32), (1, 2), (3, 5)];
    let result = contact_map_similarity(&map, &map);
    assert_eq!(result.shared, 3);
    assert_eq!(result.union, 3);
    assert!((result.jaccard - 1.0).abs() < 1e-12);
}

#[test]
fn reversed_pairs_are_the_same_contact() {
    let first = [(0u32, 1u32)];
    let second = [(1u32, 0u32)];
    let result = contact_map_similarity(&first, &second);
    assert!((result.jaccard - 1.0).abs() < 1e-12);
}

#[test]
fn disjoint_maps_share_nothing() {
    let first = [(0u32, 1u32)];
    let second = [(2u32, 3u32)];
    let result = contact_map_similarity(&first, &second);
    assert_eq!(result.shared, 0);
    assert_eq!(result.union, 2);
    assert!(result.jaccard.abs() < 1e-12);
}

#[test]
fn half_overlap_scores_one_third() {
    // {01,12} vs {12,23}: one shared, three in the union.
    let first = [(0u32, 1u32), (1, 2)];
    let second = [(1u32, 2u32), (2, 3)];
    let result = contact_map_similarity(&first, &second);
    assert_eq!(result.shared, 1);
    assert_eq!(result.union, 3);
    assert!((result.jaccard - 1.0 / 3.0).abs() < 1e-12);
}

#[test]
fn two_empty_maps_are_similar() {
    let result = contact_map_similarity(&[], &[]);
    assert!((result.jaccard - 1.0).abs() < 1e-12);
}
