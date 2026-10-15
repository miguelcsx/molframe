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

fn rdf_options() -> RadialDistributionOptions {
    RadialDistributionOptions {
        minimum_distance: 0.0,
        maximum_distance: 10.0,
        bins: 1,
        volume: 1_000.0,
        backend: SpatialBackend::BruteForce,
    }
}

/// Recovers the ideal pair population a bin was normalised by.
fn ideal_pairs(bin: &super::RadialBin, volume: f64) -> f64 {
    let shell =
        4.0 * core::f64::consts::PI * (f64::from(bin.upper).powi(3) - f64::from(bin.lower).powi(3))
            / 3.0;
    crate::numeric::u64_to_f64(bin.count) / bin.distribution / (shell / volume)
}

#[test]
fn overlapping_selections_normalize_by_unique_non_self_pairs() {
    // Sites A, B, C with left {A,B} and right {B,C}: A-B, A-C, B-C is 3 unique
    // pairs; B-B is excluded and (B,C)/(C,B) collapse.
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let Ok(bins) = radial_distribution(
        &positions,
        &selection(&[0, 1]),
        &selection(&[1, 2]),
        rdf_options(),
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid RDF");
    };
    assert_eq!(bins[0].count, 3);
    assert!((ideal_pairs(&bins[0], 1_000.0) - 3.0).abs() < 1e-9);
}

#[test]
fn identical_selections_still_use_n_choose_two() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let all = selection(&[0, 1, 2]);
    let Ok(bins) = radial_distribution(
        &positions,
        &all,
        &all,
        rdf_options(),
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid RDF");
    };
    assert_eq!(bins[0].count, 3);
    assert!((ideal_pairs(&bins[0], 1_000.0) - 3.0).abs() < 1e-9);
}

#[test]
fn centre_of_mass_overlap_never_pairs_a_group_with_itself() {
    // Three single-atom groups at 0, 1, 2 A. Left {g0, g1}, right {g1, g2}: the
    // g1-g1 pair at 0 A must not appear, leaving 3 unique pairs.
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let group = |atom: u32| CentreGroup {
        atoms: selection(&[atom]),
    };
    let left = [group(0), group(1)];
    let right = [group(1), group(2)];
    let mut options = rdf_options();
    options.bins = 10;
    let Ok(bins) = centre_of_mass_radial_distribution(
        &positions,
        &[1.0; 3],
        &left,
        &right,
        options,
        None,
        &ExecutionContext::default(),
    ) else {
        panic!("valid RDF");
    };
    assert_eq!(bins[0].count, 0, "no group may pair with itself at 0 A");
    assert_eq!(bins.iter().map(|bin| bin.count).sum::<u64>(), 3);
    let populated = bins
        .iter()
        .find(|bin| bin.count > 0)
        .unwrap_or_else(|| panic!("a populated bin"));
    // Every populated shell implies the same ideal population of 3 pairs.
    assert!((ideal_pairs(populated, 1_000.0) - 3.0).abs() < 1e-9);
}
