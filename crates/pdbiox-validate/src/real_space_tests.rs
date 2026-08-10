use super::*;
use pdbiox_core::structure::UnitCell;

fn map(values: Vec<f32>) -> DensityMap {
    DensityMap {
        dimensions: [2, 2, 1],
        starts: [0; 3],
        sampling: [2, 2, 1],
        cell: UnitCell {
            lengths: [2.0, 2.0, 1.0],
            angles: [90.0; 3],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values,
    }
}

#[test]
fn identical_maps_have_unit_correlation() {
    let observed = map(vec![0.0, 1.0, 2.0, 3.0]);
    let result = real_space_map_correlation(&observed, &observed)
        .expect("varying identical maps should correlate");
    assert!((result.coefficient - 1.0).abs() < f64::EPSILON);
    assert_eq!(result.sample_count, 4);
}

#[test]
fn mask_selects_exactly_its_voxels() {
    let observed = map(vec![0.0, 1.0, 2.0, 3.0]);
    let calculated = map(vec![0.0, 1.0, 9.0, 8.0]);
    let result = masked_real_space_correlation(&observed, &calculated, &[true, true, false, false])
        .expect("two varying masked values should correlate");
    assert!((result.coefficient - 1.0).abs() < f64::EPSILON);
    assert_eq!(result.sample_count, 2);
}

#[test]
fn sampled_correlation_omits_points_outside_both_maps() {
    let observed = map(vec![0.0, 1.0, 2.0, 3.0]);
    let calculated = map(vec![0.0, 2.0, 4.0, 6.0]);
    let result = sampled_real_space_correlation(
        &observed,
        &calculated,
        &[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [10.0, 0.0, 0.0]],
        MapBoundary::Missing,
    )
    .expect("two shared samples should correlate");
    assert_eq!(result.sample_count, 2);
    assert!((result.coefficient - 1.0).abs() < f64::EPSILON);
}

#[test]
fn constant_fields_are_reported_as_undefined() {
    let observed = map(vec![1.0; 4]);
    assert_eq!(
        real_space_map_correlation(&observed, &observed),
        Err(RealSpaceCorrelationError::ZeroVariance)
    );
}
