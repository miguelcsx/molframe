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

#[test]
fn periodic_nearest_search_uses_tree_images_across_a_triclinic_boundary() {
    let positions = [[0.1, 0.2, 0.3], [7.9, 0.2, 0.3], [4.0, 4.0, 4.0]];
    let periodic = PeriodicBox::from_cell(molframe_core::structure::UnitCell {
        lengths: [8.0, 9.0, 10.0],
        angles: [70.0, 80.0, 65.0],
    });
    let periodic = match periodic {
        Ok(periodic) => periodic,
        Err(error) => panic!("cell failed: {error}"),
    };
    let tree = match KdTree::build(&positions, &[0, 1, 2], Some(&periodic)) {
        Ok(tree) => tree,
        Err(error) => panic!("build failed: {error}"),
    };
    let nearest = match tree.k_nearest(0, 2) {
        Ok(nearest) => nearest,
        Err(error) => panic!("query failed: {error}"),
    };

    assert_eq!(nearest[0].0, 1);
    let expected = periodic.distance_squared(positions[0], positions[1]);
    assert!((nearest[0].1 - expected).abs() <= f32::EPSILON);
}

#[test]
fn periodic_image_budget_is_enforced_instead_of_falling_back() {
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
    let periodic = PeriodicBox::from_cell(molframe_core::structure::UnitCell {
        lengths: [8.0; 3],
        angles: [90.0; 3],
    });
    let periodic = match periodic {
        Ok(periodic) => periodic,
        Err(error) => panic!("cell failed: {error}"),
    };
    let tree = KdTree::build_with_options(
        &positions,
        &[0, 1],
        Some(&periodic),
        KdPeriodicOptions {
            maximum_image_count: 1,
        },
    );
    let tree = match tree {
        Ok(tree) => tree,
        Err(error) => panic!("build failed: {error}"),
    };

    assert!(matches!(
        tree.pairs(&[0], 1.0),
        Err(SpatialError::PeriodicImageLimitExceeded { .. })
    ));
}
