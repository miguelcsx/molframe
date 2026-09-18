use super::{
    CentreGroup, CoordinationOptions, RadialDistributionOptions,
    centre_of_mass_radial_distribution, coordination_numbers, radial_distribution,
};
use molframe_core::{ExecutionContext, selection::AtomSelection};
use molframe_spatial::SpatialBackend;

fn selection(indices: &[u32]) -> AtomSelection {
    AtomSelection::from_sorted(indices.to_vec())
}

#[test]
fn radial_bins_count_unique_pairs_and_normalize() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    let sites = selection(&[0, 1, 2]);
    let options = RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 4.0,
        bins: 4,
        volume: 1_000.0,
        backend: SpatialBackend::BruteForce,
    };
    let Ok(bins) = radial_distribution(
        &positions,
        &sites,
        &sites,
        options,
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid RDF");
    };

    assert_eq!(bins.iter().map(|bin| bin.count).sum::<u64>(), 3);
    assert_eq!(bins[1].count, 1);
    assert_eq!(bins[2].count, 1);
    assert_eq!(bins[3].count, 1);
    assert!(bins[1].distribution.is_finite());
}

#[test]
fn coordination_is_reported_in_left_selection_order() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
    let left = selection(&[0, 2]);
    let right = selection(&[1]);
    let Ok(counts) = coordination_numbers(
        &positions,
        &left,
        &right,
        CoordinationOptions {
            minimum_distance: 0.5,
            maximum_distance: 2.0,
            backend: SpatialBackend::BruteForce,
        },
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid coordination shell");
    };
    assert_eq!(counts, vec![1, 0]);
}

#[test]
fn radial_policy_rejects_implicit_normalization() {
    let options = RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 4.0,
        bins: 4,
        volume: 0.0,
        backend: SpatialBackend::Auto,
    };
    assert!(
        radial_distribution(
            &[],
            &selection(&[]),
            &selection(&[]),
            options,
            None,
            &ExecutionContext::default(),
        )
        .is_err()
    );
}

#[test]
fn centre_of_mass_rdf_reuses_site_distribution_over_explicit_groups() {
    let positions = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
    ];
    let groups = [
        CentreGroup {
            atoms: selection(&[0, 1]),
        },
        CentreGroup {
            atoms: selection(&[2, 3]),
        },
    ];
    let Ok(bins) = centre_of_mass_radial_distribution(
        &positions,
        &[1.0; 4],
        &groups,
        &groups,
        RadialDistributionOptions {
            minimum_distance: 0.0,
            maximum_distance: 8.0,
            bins: 4,
            volume: 1_000.0,
            backend: SpatialBackend::BruteForce,
        },
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid centre-of-mass RDF");
    };
    assert_eq!(bins.iter().map(|bin| bin.count).sum::<u64>(), 1);
    assert_eq!(bins[2].count, 1);
}
