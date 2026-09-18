use super::*;
use molframe_core::contract::MissingPolicy;
use molframe_core::structure::UnitCell;

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
fn sampled_map_coverage_records_omitted_positions() {
    let observed = map(vec![0.0, 1.0, 2.0, 3.0]);
    let calculated = map(vec![0.0, 2.0, 4.0, 6.0]);
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [10.0, 0.0, 0.0]];
    let Ok(result) = governed_sampled_real_space_correlation(
        &observed,
        &calculated,
        &positions,
        &AnalysisPolicy::default(),
    ) else {
        panic!("two varying shared samples should correlate");
    };
    assert_eq!(result.coverage.intended, 3);
    assert_eq!(result.coverage.used, 2);
    assert_eq!(result.coverage.missing, 1);
    assert_eq!(result.status, Status::Partial);
}

#[test]
fn sampled_map_obeys_fail_on_missing() {
    let observed = map(vec![0.0, 1.0, 2.0, 3.0]);
    let calculated = map(vec![0.0, 2.0, 4.0, 6.0]);
    let policy = AnalysisPolicy::default().with_missing_atoms(MissingPolicy::Fail);
    let error = governed_sampled_real_space_correlation(
        &observed,
        &calculated,
        &[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [10.0, 0.0, 0.0]],
        &policy,
    );
    assert!(matches!(error, Err(GovernedMapError::MissingSamples(1))));
}
