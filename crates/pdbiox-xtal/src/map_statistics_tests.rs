use super::*;
use pdbiox_core::structure::UnitCell;

fn map() -> DensityMap {
    DensityMap {
        dimensions: [4, 1, 1],
        starts: [0; 3],
        sampling: [4, 1, 1],
        cell: UnitCell {
            lengths: [4.0, 1.0, 1.0],
            angles: [90.0; 3],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values: vec![0.0, 1.0, 2.0, 3.0],
    }
}

#[test]
fn complete_and_masked_statistics_use_population_sigma() {
    let map = map();
    let all = map.statistics().expect("map should have statistics");
    assert!((all.mean - 1.5).abs() < f64::EPSILON);
    assert!((all.sigma - 1.25_f64.sqrt()).abs() < f64::EPSILON);
    let masked = map
        .masked_statistics(&[true, false, true, false])
        .expect("mask selects two values");
    assert!((masked.mean - 1.0).abs() < f64::EPSILON);
}

#[test]
fn histogram_includes_the_exact_upper_bound_in_last_bin() {
    let histogram = map().histogram(3, 0.0, 3.0).expect("range is valid");
    assert_eq!(histogram.counts, [1, 1, 2]);
}
