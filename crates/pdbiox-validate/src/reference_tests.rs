use super::*;
use pdbiox_core::index::{AtomIndex, ResidueIndex};

#[test]
fn scalar_assessment_retains_version_probability_and_percentile() {
    let set = ReferenceLibrary::new(
        "geometry",
        "2026.08",
        [
            ReferenceDistribution::histogram("bond_delta", vec![-1.0, 0.0, 1.0], vec![1.0, 3.0])
                .unwrap_or_else(|error| panic!("histogram failed: {error}")),
        ],
    )
    .unwrap_or_else(|error| panic!("set failed: {error}"));
    let deviation = BondDeviation {
        atom_a: AtomIndex::new(0),
        atom_b: AtomIndex::new(1),
        observed: 1.5,
        expected: 1.0,
        deviation: 0.5,
    };
    let score = assess_bond_deviation(&deviation, &set, "bond_delta")
        .unwrap_or_else(|error| panic!("assessment failed: {error}"));
    assert_eq!(score.version.as_ref(), "2026.08");
    assert!((score.probability - 0.75).abs() < f64::EPSILON);
    assert!(
        score
            .percentile
            .is_some_and(|value| (value - 0.625).abs() < f64::EPSILON)
    );
}

#[test]
fn ramachandran_grid_assessment_is_connected_to_validation_records() {
    let distribution = ReferenceDistribution::grid(
        "rama_general",
        vec![-180.0, 0.0, 180.0],
        vec![-180.0, 0.0, 180.0],
        vec![8.0, 1.0, 1.0, 0.0],
    )
    .unwrap_or_else(|error| panic!("grid failed: {error}"));
    let set = ReferenceLibrary::new("rama", "1.0", [distribution])
        .unwrap_or_else(|error| panic!("set failed: {error}"));
    let record = RamachandranRecord {
        residue: ResidueIndex::new(4),
        phi: -60.0,
        psi: -45.0,
        region: crate::RamachandranRegion::AlphaHelixRight,
        reference: ReferenceAssessment {
            set: "prior".into(),
            version: "0".into(),
            distribution: "prior".into(),
            probability: 0.0,
            percentile: None,
        },
    };
    let score = assess_ramachandran(&record, &set, "rama_general")
        .unwrap_or_else(|error| panic!("assessment failed: {error}"));
    assert!((score.probability - 0.8).abs() < f64::EPSILON);
    assert_eq!(score.distribution.as_ref(), "rama_general");
}

#[test]
fn malformed_axes_and_duplicate_names_are_rejected() {
    assert_eq!(
        ReferenceDistribution::histogram("x", vec![0.0, 0.0], vec![1.0]),
        Err(ReferenceError::Axis)
    );
    let one = ReferenceDistribution::histogram("x", vec![0.0, 1.0], vec![1.0])
        .unwrap_or_else(|error| panic!("histogram failed: {error}"));
    assert_eq!(
        ReferenceLibrary::new("set", "1", [one.clone(), one]),
        Err(ReferenceError::Duplicate)
    );
}

#[test]
fn grid_configuration_is_validated_before_observations() {
    let histogram = ReferenceDistribution::histogram("scalar", vec![0.0, 1.0], vec![1.0])
        .unwrap_or_else(|error| panic!("histogram failed: {error}"));
    let set = ReferenceLibrary::new("set", "1", [histogram])
        .unwrap_or_else(|error| panic!("library failed: {error}"));
    assert_eq!(set.validate_grid("scalar"), Err(ReferenceError::Dimension));
    assert_eq!(set.validate_grid("absent"), Err(ReferenceError::Missing));
}
