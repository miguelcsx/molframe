use super::*;
use crate::{PairQuery, for_each_pairs_within_unsorted, within};
use proptest::prelude::*;

fn context() -> molframe_core::ExecutionContext {
    molframe_core::ExecutionContext::default()
}

fn pair_indices(pairs: &[NeighborPair]) -> Vec<(u32, u32)> {
    pairs.iter().map(|pair| (pair.first, pair.second)).collect()
}

#[test]
fn within_includes_targets_and_excludes_distant_query_atoms() {
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [4.0, 0.0, 0.0]];
    let selected = within(
        &positions,
        &AtomSelection::All(3),
        &AtomSelection::from_sorted(vec![0]),
        1.0,
        SpatialBackend::CellList,
        None,
        &context(),
    );
    let selected = match selected {
        Ok(selected) => selected,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(selected.iter().collect::<Vec<_>>(), vec![0, 1]);
}

#[test]
fn every_backend_agrees_under_triclinic_periodicity() {
    let positions = [[0.1, 0.2, 0.3], [7.9, 0.2, 0.3], [4.0, 4.0, 4.0]];
    let periodic = PeriodicBox::from_cell(molframe_core::structure::UnitCell {
        lengths: [8.0, 9.0, 10.0],
        angles: [70.0, 80.0, 65.0],
    });
    let periodic = match periodic {
        Ok(periodic) => periodic,
        Err(error) => panic!("cell failed: {error}"),
    };
    let all = AtomSelection::All(3);
    let expected = match pairs_within(
        &positions,
        &all,
        &all,
        1.0,
        SpatialBackend::BruteForce,
        Some(&periodic),
        &context(),
    ) {
        Ok(pairs) => pair_indices(&pairs),
        Err(error) => panic!("query failed: {error}"),
    };
    for backend in [
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        let actual = match pairs_within(
            &positions,
            &all,
            &all,
            1.0,
            backend,
            Some(&periodic),
            &context(),
        ) {
            Ok(pairs) => pair_indices(&pairs),
            Err(error) => panic!("query failed: {error}"),
        };
        assert_eq!(actual, expected);
    }
}

#[test]
fn automatic_planning_agrees_with_brute_force_under_periodicity() {
    let positions = [[0.1, 0.2, 0.3], [7.9, 0.2, 0.3], [4.0, 4.0, 4.0]];
    let periodic = match PeriodicBox::from_cell(molframe_core::structure::UnitCell {
        lengths: [8.0, 9.0, 10.0],
        angles: [70.0, 80.0, 65.0],
    }) {
        Ok(periodic) => periodic,
        Err(error) => panic!("cell failed: {error}"),
    };
    let all = AtomSelection::All(3);
    let expected = match pairs_within_with_options(
        &positions,
        &all,
        &all,
        1.0,
        SpatialSearchOptions::with_backend(SpatialBackend::BruteForce),
        Some(&periodic),
        &context(),
    ) {
        Ok(pairs) => pair_indices(&pairs),
        Err(error) => panic!("query failed: {error}"),
    };
    let actual = match pairs_within_with_options(
        &positions,
        &all,
        &all,
        1.0,
        SpatialSearchOptions::BALANCED,
        Some(&periodic),
        &context(),
    ) {
        Ok(pairs) => pair_indices(&pairs),
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(actual, expected);
}

#[test]
fn identical_selections_keep_sorted_unique_pair_contract() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.8, 0.0, 0.0],
        [1.6, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ];
    let selection = AtomSelection::All(4);
    let expected = vec![(0, 1), (1, 2)];
    for backend in [
        SpatialBackend::BruteForce,
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        let pairs = match pairs_within(
            &positions,
            &selection,
            &selection,
            1.0,
            backend,
            None,
            &context(),
        ) {
            Ok(pairs) => pair_indices(&pairs),
            Err(error) => panic!("query failed for {backend:?}: {error}"),
        };

        assert_eq!(pairs, expected, "backend {backend:?}");
    }
}

#[test]
fn unsorted_same_selection_query_keeps_pair_set_without_order_contract() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.8, 0.0, 0.0],
        [1.6, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ];
    let selection = AtomSelection::All(4);
    let expected = pairs_within(
        &positions,
        &selection,
        &selection,
        1.0,
        SpatialBackend::CellList,
        None,
        &context(),
    )
    .unwrap_or_else(|error| panic!("sorted query failed: {error}"));
    for backend in [
        SpatialBackend::BruteForce,
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        let mut actual = pairs_within_unsorted(
            &positions,
            &selection,
            &selection,
            1.0,
            backend,
            None,
            &context(),
        )
        .unwrap_or_else(|error| panic!("unsorted query failed for {backend:?}: {error}"));
        let mut expected = expected.clone();
        actual.sort_unstable_by_key(|pair| (pair.first, pair.second));
        expected.sort_unstable_by_key(|pair| (pair.first, pair.second));
        assert_eq!(actual, expected, "backend {backend:?}");
    }
}

#[test]
fn streamed_same_selection_query_matches_materialized_pairs() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.8, 0.0, 0.0],
        [1.6, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ];
    let selection = AtomSelection::All(4);
    let expected = pairs_within_unsorted(
        &positions,
        &selection,
        &selection,
        1.0,
        SpatialBackend::CellList,
        None,
        &context(),
    )
    .unwrap_or_else(|error| panic!("materialized query failed: {error}"));
    let mut actual = Vec::new();
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions: &positions,
            left: &selection,
            right: &selection,
            cutoff: 1.0,
            options: SpatialSearchOptions::with_backend(SpatialBackend::CellList),
            periodic: None,
            context: &context(),
        },
        |pair| actual.push(pair),
    )
    .unwrap_or_else(|error| panic!("streamed query failed: {error}"));
    actual.sort_unstable_by_key(|pair| (pair.first, pair.second));
    let mut expected = expected;
    expected.sort_unstable_by_key(|pair| (pair.first, pair.second));
    assert_eq!(actual, expected);
}

#[test]
fn streamed_overlapping_cross_selections_keep_the_unique_pair_contract() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.6, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [1.8, 0.0, 0.0],
        [2.4, 0.0, 0.0],
    ];
    let left = AtomSelection::from_sorted(vec![0, 1, 2, 4]);
    let right = AtomSelection::from_sorted(vec![1, 2, 3, 4]);
    for backend in [SpatialBackend::BruteForce, SpatialBackend::CellList] {
        let mut expected =
            pairs_within_unsorted(&positions, &left, &right, 1.3, backend, None, &context())
                .unwrap_or_else(|error| {
                    panic!("materialized query failed for {backend:?}: {error}")
                });
        let stream = || {
            let mut pairs = Vec::new();
            for_each_pairs_within_unsorted(
                &PairQuery {
                    positions: &positions,
                    left: &left,
                    right: &right,
                    cutoff: 1.3,
                    options: SpatialSearchOptions::with_backend(backend),
                    periodic: None,
                    context: &context(),
                },
                |pair| pairs.push(pair),
            )
            .unwrap_or_else(|error| panic!("streamed query failed for {backend:?}: {error}"));
            pairs
        };
        let mut actual = stream();
        assert_eq!(actual, stream(), "stream order changed for {backend:?}");
        expected.sort_unstable_by_key(|pair| (pair.first, pair.second));
        actual.sort_unstable_by_key(|pair| (pair.first, pair.second));
        assert_eq!(actual, expected, "{backend:?}");
    }
}

proptest! {
    #[test]
    fn all_backends_return_the_same_pairs(
        positions in prop::collection::vec(
            (-20.0_f32..20.0, -20.0_f32..20.0, -20.0_f32..20.0),
            1..80,
        ),
        cutoff in 0.0_f32..8.0,
    ) {
        let positions: Vec<[f32; 3]> = positions
            .into_iter()
            .map(|(x, y, z)| [x, y, z])
            .collect();
        let Ok(atom_count) = u32::try_from(positions.len()) else {
            panic!("test position count fits u32");
        };
        let all = AtomSelection::All(atom_count);
        let expected = pairs_within(
            &positions,
            &all,
            &all,
            cutoff,
            SpatialBackend::BruteForce,
            None,
            &context(),
        );
        let expected = match expected {
            Ok(pairs) => pair_indices(&pairs),
            Err(error) => return Err(TestCaseError::fail(error.to_string())),
        };
        for backend in [
            SpatialBackend::CellList,
            SpatialBackend::KdTree,
            SpatialBackend::NeighborList,
        ] {
            let actual = pairs_within(
                &positions, &all, &all, cutoff, backend, None, &context(),
            );
            let actual = match actual {
                Ok(pairs) => pair_indices(&pairs),
                Err(error) => return Err(TestCaseError::fail(error.to_string())),
            };
            prop_assert_eq!(actual, expected.clone());
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]
    #[test]
    fn native_periodic_indices_match_minimum_image_brute_force(
        positions in prop::collection::vec(
            (-12.0_f32..12.0, -12.0_f32..12.0, -12.0_f32..12.0),
            1..35,
        ),
        cutoff in 0.0_f32..3.0,
        triclinic in any::<bool>(),
    ) {
        let positions: Vec<[f32; 3]> = positions
            .into_iter()
            .map(|(x, y, z)| [x, y, z])
            .collect();
        let angles = if triclinic {
            [70.0, 80.0, 65.0]
        } else {
            [90.0; 3]
        };
        let periodic = PeriodicBox::from_cell(molframe_core::structure::UnitCell {
            lengths: [8.0, 9.0, 10.0],
            angles,
        });
        let periodic = match periodic {
            Ok(periodic) => periodic,
            Err(error) => return Err(TestCaseError::fail(error.to_string())),
        };
        let Ok(atom_count) = u32::try_from(positions.len()) else {
            panic!("test position count fits u32");
        };
        let all = AtomSelection::All(atom_count);
        let expected = pairs_within(
            &positions,
            &all,
            &all,
            cutoff,
            SpatialBackend::BruteForce,
            Some(&periodic),
            &context(),
        );
        let expected = match expected {
            Ok(pairs) => pair_indices(&pairs),
            Err(error) => return Err(TestCaseError::fail(error.to_string())),
        };

        for backend in [SpatialBackend::CellList, SpatialBackend::KdTree] {
            let actual = pairs_within(
                &positions,
                &all,
                &all,
                cutoff,
                backend,
                Some(&periodic),
                &context(),
            );
            let actual = match actual {
                Ok(pairs) => pair_indices(&pairs),
                Err(error) => return Err(TestCaseError::fail(error.to_string())),
            };
            prop_assert_eq!(actual, expected.clone());
        }
    }
}

// A dense cubic lattice forces the cell-list backend and many non-empty cells,
// so the parallel partition spans multiple cells per worker.
fn lattice(side: u8) -> Vec<[f32; 3]> {
    let mut positions = Vec::new();
    for x in 0..side {
        for y in 0..side {
            for z in 0..side {
                positions.push([f32::from(x), f32::from(y), f32::from(z)]);
            }
        }
    }
    positions
}

#[test]
fn parallel_pairs_match_serial_for_every_worker_count() {
    let positions = lattice(16); // 4096 atoms -> cell-list backend
    let all = AtomSelection::All(4096);
    let options = SpatialSearchOptions::with_backend(SpatialBackend::CellList);
    let serial = pairs_within_with_options(&positions, &all, &all, 1.5, options, None, &context())
        .expect("serial cell-list search is valid");
    for workers in [1, 2, 3, 4, 8, 16] {
        let worker_context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("worker context is valid");
        let parallel =
            pairs_within_with_options(&positions, &all, &all, 1.5, options, None, &worker_context)
                .expect("parallel cell-list search is valid");
        assert_eq!(
            pair_indices(&serial),
            pair_indices(&parallel),
            "worker count {workers} changed the pair set"
        );
        // Distances must match too, not only the index pairs.
        assert_eq!(
            serial
                .iter()
                .map(|p| p.distance_squared.to_bits())
                .collect::<Vec<_>>(),
            parallel
                .iter()
                .map(|p| p.distance_squared.to_bits())
                .collect::<Vec<_>>(),
            "worker count {workers} changed a squared distance"
        );
    }
}

#[test]
fn sparse_index_validation_borrows_the_existing_allocation() {
    let selection = AtomSelection::Sparse(vec![1, 5, 9]);
    let indices = super::checked_indices(&selection, 10).expect("valid indices");
    let AtomSelection::Sparse(original) = &selection else {
        panic!("sparse input")
    };
    assert!(matches!(indices, std::borrow::Cow::Borrowed(_)));
    assert_eq!(indices.as_ptr(), original.as_ptr());
    assert_eq!(
        super::checked_indices(&selection, 9),
        Err(SpatialError::AtomOutOfBounds(9))
    );
}
