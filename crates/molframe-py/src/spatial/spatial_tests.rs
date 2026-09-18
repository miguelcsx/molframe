use super::*;

#[test]
fn neighbor_pair_is_a_constructible_canonical_scalar() {
    let pair = match PyNeighborPair::new(2, 7, 3.5) {
        Ok(pair) => pair,
        Err(error) => panic!("ascending finite pair should construct: {error}"),
    };
    assert_eq!(pair.first, 2);
    assert_eq!(pair.second, 7);
    assert!((pair.distance - 3.5).abs() < f32::EPSILON);
    assert_eq!(
        pair.__repr__(),
        "NeighborPair(first=2, second=7, distance=3.5)"
    );
}

#[test]
fn neighbor_pair_rejects_noncanonical_or_invalid_values() {
    assert!(PyNeighborPair::new(7, 2, 3.5).is_err());
    assert!(PyNeighborPair::new(2, 2, 3.5).is_err());
    assert!(PyNeighborPair::new(2, 7, -0.1).is_err());
    assert!(PyNeighborPair::new(2, 7, f32::NAN).is_err());
}
