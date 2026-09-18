use super::{Linkage, agglomerative_clustering, dbscan_clustering, medoid};
use crate::EnsembleDistanceMatrix;

fn distances() -> EnsembleDistanceMatrix {
    EnsembleDistanceMatrix {
        size: 4,
        values: vec![
            0.0, 0.2, 5.0, 5.1, 0.2, 0.0, 4.9, 5.0, 5.0, 4.9, 0.0, 0.1, 5.1, 5.0, 0.1, 0.0,
        ]
        .into_boxed_slice(),
    }
}

#[test]
fn agglomerative_and_dbscan_find_the_two_dense_groups() {
    let matrix = distances();
    let Ok(hierarchical) = agglomerative_clustering(&matrix, 2, Linkage::Average) else {
        panic!("valid clustering");
    };
    let Ok(density) = dbscan_clustering(&matrix, 0.25, 2) else {
        panic!("valid density clustering");
    };
    assert_eq!(hierarchical.members, vec![vec![0, 1], vec![2, 3]]);
    assert_eq!(density.members, hierarchical.members);
}

#[test]
fn medoid_ties_choose_the_lowest_observation() {
    let matrix = distances();
    assert_eq!(medoid(&matrix, &[2, 3]), Ok(2));
}
