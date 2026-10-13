use super::PairColumns;
use crate::NeighborPair;

#[test]
fn columns_hold_indices_and_euclidean_distances() {
    let pairs = [NeighborPair::new(3, 1, 4.0), NeighborPair::new(0, 2, 9.0)];
    let columns = PairColumns::from_pairs(&pairs);
    assert_eq!(columns.first, [1, 0]);
    assert_eq!(columns.second, [3, 2]);
    assert_eq!(columns.distance, [2.0, 3.0]);
}
