use super::{
    AnisotropicNetworkModel, AnmError, AnmOptions,
    anisotropic_network_model as anisotropic_network_model_native,
};
use crate::numeric::f64_to_f32;
use molframe_core::{ExecutionContext, selection::AtomSelection};
use molframe_spatial::{PeriodicBox, SpatialBackend};

fn anisotropic_network_model(
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    options: AnmOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<AnisotropicNetworkModel, AnmError> {
    anisotropic_network_model_native(
        positions,
        sites,
        options,
        periodic,
        &ExecutionContext::default(),
    )
}

fn options(modes: usize) -> AnmOptions {
    AnmOptions {
        contact_distance: 7.5,
        mode_count: modes,
        zero_mode_tolerance: 1e-8,
        memory_limit_bytes: 64 * 1_024 * 1_024,
        backend: SpatialBackend::BruteForce,
    }
}

fn sites(count: u32) -> AtomSelection {
    AtomSelection::from_sorted((0..count).collect())
}

/// Alpha-carbon positions along an ideal alpha helix.
///
/// A straight line of sites is deliberately *not* used: springs along a
/// collinear chain only ever pull along one axis, which leaves the other two
/// completely unconstrained and makes the network degenerate. A helix is what
/// a real backbone looks like and what the cutoff was chosen for.
fn helix(count: usize) -> Vec<[f32; 3]> {
    const RADIUS: f64 = 2.3;
    const RISE: f64 = 1.5;
    const TURN: f64 = 1.745_329_25;
    (0..count)
        .map(|index| {
            let Ok(index) = u16::try_from(index) else {
                panic!("helix fixture fits u16")
            };
            let step = f64::from(index);
            let angle = step * TURN;
            [
                f64_to_f32(RADIUS * angle.cos()),
                f64_to_f32(RADIUS * angle.sin()),
                f64_to_f32(RISE * step),
            ]
        })
        .collect()
}

#[test]
fn a_helix_reports_exactly_the_six_rigid_directions() {
    let positions = helix(6);
    let Ok(model) = anisotropic_network_model(&positions, &sites(6), options(3), None) else {
        panic!("valid connected network")
    };
    // Three translations and three rotations, and nothing else: a helix has
    // extent about every axis, so no further direction is free of restoring
    // force.
    assert_eq!(model.zero_modes, 6);
    assert_eq!(model.eigenvalues.len(), 3);
    assert_eq!(model.modes.len(), 3);
    assert_eq!(model.modes[0].len(), 6);
}

#[test]
fn a_collinear_chain_is_degenerate_beyond_its_rigid_directions() {
    // Springs along a straight line only ever resist motion along that line,
    // so both transverse axes are free at every site. The report counts those
    // directions honestly rather than pretending the network constrains them.
    let positions: Vec<[f32; 3]> = (0_u16..6)
        .map(|index| [f32::from(index) * 3.8, 0.0, 0.0])
        .collect();
    let Ok(model) = anisotropic_network_model(&positions, &sites(6), options(3), None) else {
        panic!("valid connected network")
    };
    assert_eq!(model.zero_modes, 13);
}

#[test]
fn eigenvalues_come_back_in_increasing_order_and_are_positive() {
    let positions = helix(8);
    let Ok(model) = anisotropic_network_model(&positions, &sites(8), options(4), None) else {
        panic!("valid connected network")
    };
    for value in &model.eigenvalues {
        assert!(*value > 0.0, "a non-zero mode has positive stiffness");
    }
    for pair in model.eigenvalues.windows(2) {
        assert!(pair[0] <= pair[1], "eigenvalues are sorted: {pair:?}");
    }
}

#[test]
fn modes_are_unit_length_and_free_of_net_translation() {
    let positions = helix(10);
    let Ok(model) = anisotropic_network_model(&positions, &sites(10), options(4), None) else {
        panic!("valid connected network")
    };
    for mode in &model.modes {
        let norm: f64 = mode
            .iter()
            .flatten()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        assert!((norm - 1.0).abs() < 1e-6, "mode norm was {norm}");
        for axis in 0..3 {
            let drift: f64 = mode.iter().map(|displacement| displacement[axis]).sum();
            assert!(drift.abs() < 1e-6, "mode drifts on axis {axis}: {drift}");
        }
    }
}

#[test]
fn the_same_input_always_returns_the_same_modes() {
    let positions = helix(12);
    let Ok(first) = anisotropic_network_model(&positions, &sites(12), options(3), None) else {
        panic!("valid connected network")
    };
    let Ok(second) = anisotropic_network_model(&positions, &sites(12), options(3), None) else {
        panic!("valid connected network")
    };
    assert_eq!(first, second);
}

#[test]
fn projecting_a_mode_onto_itself_recovers_its_amplitude() {
    let positions = helix(9);
    let Ok(model) = anisotropic_network_model(&positions, &sites(9), options(3), None) else {
        panic!("valid connected network")
    };
    let scaled: Vec<[f64; 3]> = model.modes[0]
        .iter()
        .map(|displacement| displacement.map(|value| value * 2.5))
        .collect();
    let overlaps = model.project(&scaled);
    assert!(
        (overlaps[0] - 2.5).abs() < 1e-6,
        "self-overlap was {}",
        overlaps[0]
    );
    for overlap in &overlaps[1..] {
        assert!(overlap.abs() < 1e-6, "modes are orthogonal: {overlap}");
    }
}

#[test]
fn displacing_along_a_mode_moves_every_site_by_its_own_vector() {
    let positions = helix(7);
    let Ok(model) = anisotropic_network_model(&positions, &sites(7), options(2), None) else {
        panic!("valid connected network")
    };
    let mut moved = vec![[0.0_f32; 3]; positions.len()];
    model.displace(&positions, &[3.0], &mut moved);
    for ((original, shifted), displacement) in positions.iter().zip(&moved).zip(&model.modes[0]) {
        for axis in 0..3 {
            let expected = f64::from(original[axis]) + 3.0 * displacement[axis];
            assert!(
                (f64::from(shifted[axis]) - expected).abs() < 1e-4,
                "axis {axis} moved to {} not {expected}",
                shifted[axis]
            );
        }
    }
}

#[test]
fn fluctuations_are_largest_where_a_chain_is_least_constrained() {
    let positions = helix(11);
    let Ok(model) = anisotropic_network_model(&positions, &sites(11), options(6), None) else {
        panic!("valid connected network")
    };
    let fluctuations = model.fluctuations();
    assert_eq!(fluctuations.len(), positions.len());
    let middle = fluctuations.len() / 2;
    assert!(
        fluctuations[0] > fluctuations[middle],
        "chain ends whip more than its middle: {fluctuations:?}"
    );
}

#[test]
fn a_protein_sized_network_solves_and_reports_six_rigid_directions() {
    // Two hundred residues is an ordinary single domain, and the size the
    // direct solve is meant to be comfortable at.
    let positions = helix(200);
    let Ok(model) = anisotropic_network_model(&positions, &sites(200), options(10), None) else {
        panic!("valid connected network")
    };
    assert_eq!(model.zero_modes, 6);
    assert_eq!(model.eigenvalues.len(), 10);
    assert_eq!(model.modes[0].len(), 200);
    for pair in model.eigenvalues.windows(2) {
        assert!(pair[0] <= pair[1], "eigenvalues are sorted: {pair:?}");
    }
}

#[test]
fn a_network_too_large_for_the_ceiling_is_refused_rather_than_attempted() {
    let positions = helix(200);
    let mut constrained = options(1);
    // Far more than the site table needs, far less than a 600x600 dense
    // matrix: the refusal must name the matrix, not the selection.
    constrained.memory_limit_bytes = 64 * 1_024;
    assert!(matches!(
        anisotropic_network_model(&positions, &sites(200), constrained, None),
        Err(AnmError::MemoryLimit { required, limit })
            if required > limit && limit == 64 * 1_024
    ));
}

#[test]
fn a_caller_ceiling_is_refused_before_anything_large_is_reserved() {
    let positions = helix(4);
    let mut constrained = options(1);
    let exact = 4 * size_of::<u32>();
    constrained.memory_limit_bytes = exact - 1;
    assert!(matches!(
        anisotropic_network_model(&positions, &sites(4), constrained, None),
        Err(AnmError::MemoryLimit { required, limit })
            if required == exact && limit == exact - 1
    ));
}

#[test]
fn invalid_controls_and_short_selections_are_refused() {
    let positions = helix(4);
    let mut invalid = options(1);
    invalid.contact_distance = 0.0;
    assert!(matches!(
        anisotropic_network_model(&positions, &sites(4), invalid, None),
        Err(AnmError::InvalidOptions)
    ));
    assert!(matches!(
        anisotropic_network_model(&positions, &sites(1), options(1), None),
        Err(AnmError::TooFewSites)
    ));
}

#[test]
fn more_modes_than_the_network_has_are_refused() {
    let positions = helix(3);
    // Nine degrees of freedom less six rigid ones leaves three.
    assert!(matches!(
        anisotropic_network_model(&positions, &sites(3), options(4), None),
        Err(AnmError::InsufficientModes { requested: 4, .. })
    ));
}
