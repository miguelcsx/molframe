use super::{LeafletOptions, identify_leaflets};
use pdbiox_core::{ExecutionContext, selection::AtomSelection};
use pdbiox_spatial::SpatialBackend;

#[test]
fn separated_site_graphs_form_two_leaflets() {
    let positions = [
        [0.0, 0.0, -2.0],
        [1.0, 0.0, -2.0],
        [0.0, 0.0, 2.0],
        [1.0, 0.0, 2.0],
    ];
    let sites = AtomSelection::from_sorted(vec![0, 1, 2, 3]);
    let Ok(leaflets) = identify_leaflets(
        &positions,
        &sites,
        LeafletOptions {
            connection_distance: 1.5,
            backend: SpatialBackend::BruteForce,
        },
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid leaflet graph");
    };
    assert_eq!(leaflets.len(), 2);
    assert_eq!(leaflets[0].sites, vec![0, 1]);
    assert_eq!(leaflets[1].sites, vec![2, 3]);
}

#[test]
fn caller_selection_controls_chemical_membership() {
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let sites = AtomSelection::from_sorted(vec![0, 2]);
    let Ok(leaflets) = identify_leaflets(
        &positions,
        &sites,
        LeafletOptions {
            connection_distance: 2.0,
            backend: SpatialBackend::BruteForce,
        },
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid selection");
    };
    assert_eq!(leaflets[0].sites, vec![0, 2]);
}
