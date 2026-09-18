use super::polymer_statistics;

#[test]
fn bent_path_reports_contour_end_to_end_and_persistence() {
    let path = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.866_025_4, 0.0]];
    let Ok(statistics) = polymer_statistics(&path) else {
        panic!("valid polymer path");
    };
    assert!((statistics.contour_length - 2.0).abs() < 1e-6);
    let Some(correlation) = statistics.adjacent_tangent_correlation else {
        panic!("bent path should have a tangent correlation");
    };
    assert!((correlation - 0.5).abs() < 1e-6);
    assert!(statistics.persistence_length.is_some());
}

#[test]
fn straight_path_does_not_fabricate_an_infinite_persistence_length() {
    let Ok(statistics) = polymer_statistics(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]) else {
        panic!("valid straight path");
    };
    assert_eq!(statistics.adjacent_tangent_correlation, Some(1.0));
    assert_eq!(statistics.persistence_length, None);
}

#[test]
fn polymer_order_must_be_explicit_and_non_degenerate() {
    assert!(polymer_statistics(&[[0.0; 3], [0.0; 3]]).is_err());
}
