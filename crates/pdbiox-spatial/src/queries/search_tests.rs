use super::*;
use proptest::prelude::*;

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
    let periodic = PeriodicBox::from_cell(pdbiox_core::structure::UnitCell {
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
    ) {
        Ok(pairs) => pair_indices(&pairs),
        Err(error) => panic!("query failed: {error}"),
    };
    for backend in [
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::NeighborList,
    ] {
        let actual = match pairs_within(&positions, &all, &all, 1.0, backend, Some(&periodic)) {
            Ok(pairs) => pair_indices(&pairs),
            Err(error) => panic!("query failed: {error}"),
        };
        assert_eq!(actual, expected);
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
            let actual = pairs_within(&positions, &all, &all, cutoff, backend, None);
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
        let periodic = PeriodicBox::from_cell(pdbiox_core::structure::UnitCell {
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
            );
            let actual = match actual {
                Ok(pairs) => pair_indices(&pairs),
                Err(error) => return Err(TestCaseError::fail(error.to_string())),
            };
            prop_assert_eq!(actual, expected.clone());
        }
    }
}
