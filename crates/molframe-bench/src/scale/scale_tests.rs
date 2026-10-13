//! Determinism and arithmetic of the scale-sweep infrastructure.

use super::*;

#[test]
fn exponent_of_an_exact_power_law_is_recovered() {
    let samples: Vec<(f64, f64)> = [1.0e3, 1.0e4, 1.0e5, 1.0e6]
        .iter()
        .map(|&n: &f64| (n, 3.0e-9 * n.powf(1.5)))
        .collect();
    let fit = fit_exponent(&samples).unwrap();
    assert!((fit.exponent - 1.5).abs() < 1e-9);
    assert!(fit.r_squared > 0.999_999);
}

#[test]
fn exponent_needs_two_positive_samples() {
    assert!(fit_exponent(&[(10.0, 1.0)]).is_none());
    assert!(fit_exponent(&[(10.0, 1.0), (10.0, 2.0)]).is_none());
    assert!(fit_exponent(&[(10.0, 0.0), (20.0, 1.0)]).is_none());
}

#[test]
fn generators_are_pure_functions_of_the_seed() {
    assert_eq!(chain(500, 7), chain(500, 7));
    assert_ne!(chain(500, 7), chain(500, 8));
    assert_eq!(
        uniform_cloud(500, 3).positions,
        uniform_cloud(500, 3).positions
    );
    assert_eq!(
        ensemble(20, 5, 0.5, 1).positions,
        ensemble(20, 5, 0.5, 1).positions
    );
}

#[test]
fn chain_steps_are_one_bond_long() {
    let points = chain(200, 11);
    for pair in points.windows(2) {
        let d: f32 = (0..3)
            .map(|a| (pair[1][a] - pair[0][a]).powi(2))
            .sum::<f32>()
            .sqrt();
        assert!((d - 3.8).abs() < 1e-3, "{d}");
    }
}

#[test]
fn rotation_preserves_internal_distances() {
    let points = chain(100, 5);
    let moved = rotated(&points, 9);
    let d = |p: &[[f32; 3]], i: usize, j: usize| -> f32 {
        (0..3)
            .map(|a| (p[i][a] - p[j][a]).powi(2))
            .sum::<f32>()
            .sqrt()
    };
    assert!((d(&points, 3, 77) - d(&moved, 3, 77)).abs() < 1e-2);
}

#[test]
fn uniform_cloud_density_is_protein_like() {
    let cloud = uniform_cloud(8_000, 2);
    let density = 8_000.0 / cloud.edge.powi(3);
    assert!((density - 0.08).abs() < 1e-3);
}

#[test]
fn checksum_is_bit_and_order_sensitive() {
    let mut a = Checksum::default();
    a.reals(&[1.0, 2.0]);
    let mut b = Checksum::default();
    b.reals(&[2.0, 1.0]);
    let mut c = Checksum::default();
    c.reals(&[1.0, f64::from_bits(2.0_f64.to_bits() + 1)]);
    assert_ne!(a, b);
    assert_ne!(a, c);
}

#[test]
fn tiers_round_trip_through_their_labels() {
    for tier in Tier::ALL {
        assert_eq!(Tier::parse(tier.label()).unwrap(), tier);
    }
    assert!(Tier::parse("gigantic").is_err());
    assert!(Tier::Tiny < Tier::Huge);
}
