use super::solver::modes_from_matrix;
use super::{
    GaussianNetworkModel, GnmError, GnmOptions,
    gaussian_network_model as gaussian_network_model_native,
};
use molframe_core::{ExecutionContext, selection::AtomSelection};
use molframe_spatial::{PeriodicBox, SpatialBackend, pairs_within};
use nalgebra::DMatrix;
use std::collections::BTreeMap;

fn options(modes: usize) -> GnmOptions {
    GnmOptions {
        contact_distance: 1.1,
        mode_count: modes,
        zero_mode_tolerance: 1e-10,
        memory_limit_bytes: 64 * 1_024 * 1_024,
        backend: SpatialBackend::BruteForce,
        reduction: molframe_core::parallel::ReductionPolicy::Deterministic,
    }
}

fn gaussian_network_model(
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    options: GnmOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<GaussianNetworkModel, GnmError> {
    gaussian_network_model_native(
        positions,
        sites,
        options,
        periodic,
        &ExecutionContext::default(),
    )
}

#[test]
fn connected_three_site_network_has_one_zero_mode() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let sites = AtomSelection::from_sorted(vec![0, 1, 2]);
    let Ok(model) = gaussian_network_model(&positions, &sites, options(2), None) else {
        panic!("valid connected network");
    };
    assert_eq!(model.zero_modes, 1);
    assert_eq!(model.eigenvalues.len(), 2);
    assert!(model.eigenvalues[0] > 0.0);
}

#[test]
fn caller_controls_memory_and_zero_mode_policy() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let sites = AtomSelection::from_sorted(vec![0, 1]);
    let mut constrained = options(1);
    let exact = 2 * size_of::<u32>();
    constrained.memory_limit_bytes = exact - 1;
    assert!(matches!(
        gaussian_network_model(&positions, &sites, constrained, None),
        Err(GnmError::MemoryLimit { required, limit })
            if required == exact && limit == exact - 1
    ));
    constrained.memory_limit_bytes = 4_096;
    assert!(gaussian_network_model(&positions, &sites, constrained, None).is_ok());
}

#[test]
fn sparse_modes_match_the_dense_reference() {
    let positions: Vec<[f32; 3]> = (0_u16..96)
        .map(|index| [f32::from(index), 0.0, 0.0])
        .collect();
    let sites = AtomSelection::All(96);
    let current = options(4);
    let Ok(actual) = gaussian_network_model(&positions, &sites, current, None) else {
        panic!("sparse path network must solve");
    };
    let Ok(expected) = materialized_reference(&positions, &sites, current) else {
        panic!("dense path reference must solve");
    };
    assert_models_close(&actual, &expected, 1e-8);
}

#[test]
fn sparse_solver_applies_the_explicit_numerical_zero_tolerance() {
    let positions: Vec<[f32; 3]> = (0_u16..96)
        .map(|index| [f32::from(index), 0.0, 0.0])
        .collect();
    let sites = AtomSelection::All(96);
    let mut current = options(2);
    current.zero_mode_tolerance = 0.002;
    let Ok(actual) = gaussian_network_model(&positions, &sites, current, None) else {
        panic!("sparse path with a numerical zero must solve");
    };
    let Ok(expected) = materialized_reference(&positions, &sites, current) else {
        panic!("dense numerical-zero reference must solve");
    };
    assert_eq!(actual.zero_modes, 2);
    assert_models_close(&actual, &expected, 1e-8);
}

#[test]
fn sparse_disconnected_network_counts_every_component_zero_mode() {
    let first = (0_u16..41).map(|index| [f32::from(index), 0.0, 0.0]);
    let second = (0_u16..43).map(|index| [100.0 + f32::from(index), 0.0, 0.0]);
    let positions: Vec<[f32; 3]> = first.chain(second).collect();
    let sites = AtomSelection::All(84);
    let current = options(4);
    let Ok(actual) = gaussian_network_model(&positions, &sites, current, None) else {
        panic!("disconnected sparse network must solve");
    };
    let Ok(expected) = materialized_reference(&positions, &sites, current) else {
        panic!("disconnected dense reference must solve");
    };
    assert_eq!(actual.zero_modes, 2);
    assert_models_close(&actual, &expected, 1e-7);
}

#[test]
fn sparse_solver_recovers_repeated_modes_across_disconnected_components() {
    let mut positions = Vec::with_capacity(66);
    for component in 0_u16..33 {
        let origin = 3.0 * f32::from(component);
        positions.push([origin, 0.0, 0.0]);
        positions.push([origin + 1.0, 0.0, 0.0]);
    }
    let sites = AtomSelection::All(66);
    let current = options(3);
    let Ok(model) = gaussian_network_model(&positions, &sites, current, None) else {
        panic!("repeated sparse modes must converge through deflation");
    };
    assert_eq!(model.zero_modes, 33);
    for eigenvalue in &model.eigenvalues {
        assert!((eigenvalue - 2.0).abs() <= 1e-10);
    }
    for left in 0..model.modes.len() {
        for right in 0..left {
            let overlap: f64 = model.modes[left]
                .iter()
                .zip(&model.modes[right])
                .map(|(left, right)| left * right)
                .sum();
            assert!(overlap.abs() <= 1e-10);
        }
    }
}

#[test]
fn streamed_contacts_match_the_materialized_reference() {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.75, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.25, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ];
    let sites = AtomSelection::All(5);
    for backend in [SpatialBackend::BruteForce, SpatialBackend::CellList] {
        let mut current = options(4);
        current.backend = backend;
        let expected = materialized_reference(&positions, &sites, current);
        let actual = gaussian_network_model(&positions, &sites, current, None);
        match (actual, expected) {
            (Ok(actual), Ok(expected)) => assert_eq!(actual, expected, "{backend:?}"),
            (actual, expected) => {
                panic!("GNM equivalence failed for {backend:?}: {actual:?} {expected:?}")
            }
        }
    }
}

#[test]
fn ten_thousand_site_grid_solves_inside_the_extreme_memory_ceiling() {
    let mut positions = Vec::with_capacity(10_000);
    for z in 0_u16..20 {
        for y in 0_u16..20 {
            for x in 0_u16..25 {
                positions.push([f32::from(x), f32::from(y), f32::from(z)]);
            }
        }
    }
    let sites = AtomSelection::All(10_000);
    let mut current = options(1);
    current.contact_distance = 1.01;
    current.backend = SpatialBackend::CellList;
    current.memory_limit_bytes = 500_000_000;
    let Ok(model) = gaussian_network_model(&positions, &sites, current, None) else {
        panic!("10k sparse GNM must solve rather than reject a dense workspace");
    };
    assert_eq!(model.zero_modes, 1);
    assert_eq!(model.eigenvalues.len(), 1);
    assert_eq!(model.modes[0].len(), 10_000);
    assert!(model.eigenvalues[0] > current.zero_mode_tolerance);
}

fn assert_models_close(
    actual: &GaussianNetworkModel,
    expected: &GaussianNetworkModel,
    tolerance: f64,
) {
    assert_eq!(actual.sites, expected.sites);
    assert_eq!(actual.zero_modes, expected.zero_modes);
    assert_eq!(actual.eigenvalues.len(), expected.eigenvalues.len());
    for (actual_value, expected_value) in actual.eigenvalues.iter().zip(&expected.eigenvalues) {
        assert!(
            (actual_value - expected_value).abs() <= tolerance,
            "eigenvalue mismatch: {actual_value} != {expected_value}"
        );
    }
    for (actual_mode, expected_mode) in actual.modes.iter().zip(&expected.modes) {
        let overlap: f64 = actual_mode
            .iter()
            .zip(expected_mode)
            .map(|(actual, expected)| actual * expected)
            .sum();
        assert!(
            overlap.abs() >= 1.0 - tolerance,
            "mode overlap {overlap} is below tolerance"
        );
    }
}

fn materialized_reference(
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    options: GnmOptions,
) -> Result<GaussianNetworkModel, GnmError> {
    let selected: Vec<u32> = sites.into_iter().collect();
    let contacts = pairs_within(
        positions,
        sites,
        sites,
        options.contact_distance,
        options.backend,
        None,
        &ExecutionContext::default(),
    )?;
    let lookup: BTreeMap<u32, usize> = selected
        .iter()
        .copied()
        .enumerate()
        .map(|(local, atom)| (atom, local))
        .collect();
    let mut kirchhoff = DMatrix::zeros(selected.len(), selected.len());
    for pair in contacts {
        let Some(&left) = lookup.get(&pair.first) else {
            continue;
        };
        let Some(&right) = lookup.get(&pair.second) else {
            continue;
        };
        kirchhoff[(left, right)] = -1.0;
        kirchhoff[(right, left)] = -1.0;
        kirchhoff[(left, left)] += 1.0;
        kirchhoff[(right, right)] += 1.0;
    }
    modes_from_matrix(kirchhoff, selected, options)
}

/// A straight chain of unit-spaced sites, connected under the test cutoff.
fn chain(count: u32) -> (Vec<[f32; 3]>, AtomSelection) {
    let positions = (0..count)
        .map(|index| {
            [
                f32::from(u16::try_from(index).expect("small chain")),
                0.0,
                0.0,
            ]
        })
        .collect();
    (positions, AtomSelection::All(count))
}

#[test]
fn the_deterministic_policy_gives_identical_eigenvalues_on_every_solve() {
    let (positions, sites) = chain(24);
    let Ok(first) = gaussian_network_model(&positions, &sites, options(4), None) else {
        panic!("valid chain network");
    };
    let Ok(second) = gaussian_network_model(&positions, &sites, options(4), None) else {
        panic!("valid chain network");
    };

    assert_eq!(
        first.eigenvalues, second.eigenvalues,
        "the deterministic policy must be bit-identical across solves"
    );
}

#[test]
fn the_fast_policy_agrees_with_the_deterministic_one_within_solver_tolerance() {
    // The relaxed policy lets the dense solver reorder its reductions, so the
    // eigenvalues agree to the solver's tolerance rather than exactly.
    let (positions, sites) = chain(24);
    let Ok(exact) = gaussian_network_model(&positions, &sites, options(4), None) else {
        panic!("valid chain network");
    };

    let mut relaxed_options = options(4);
    relaxed_options.reduction = molframe_core::parallel::ReductionPolicy::Fast;
    let Ok(relaxed) = gaussian_network_model(&positions, &sites, relaxed_options, None) else {
        panic!("valid chain network");
    };

    assert_eq!(relaxed.eigenvalues.len(), exact.eigenvalues.len());
    for (fast, reference) in relaxed.eigenvalues.iter().zip(&exact.eigenvalues) {
        assert!(
            (fast - reference).abs() <= 1.0e-6 * reference.abs().max(1.0),
            "relaxed eigenvalue {fast} differs from {reference} beyond tolerance"
        );
    }
}
