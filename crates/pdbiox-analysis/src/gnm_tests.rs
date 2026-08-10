use super::{GnmOptions, gaussian_network_model};
use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::SpatialBackend;

fn options(modes: usize) -> GnmOptions {
    GnmOptions {
        contact_distance: 1.1,
        mode_count: modes,
        zero_mode_tolerance: 1e-10,
        memory_limit_bytes: 1_024,
        backend: SpatialBackend::BruteForce,
    }
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
    constrained.memory_limit_bytes = 1;
    assert!(gaussian_network_model(&positions, &sites, constrained, None).is_err());
}
