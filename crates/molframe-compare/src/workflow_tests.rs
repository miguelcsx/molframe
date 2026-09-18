use super::*;

#[test]
fn four_stages_retain_mapping_transform_raw_values_and_verdict() {
    let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let model = [[5.0, 2.0, -1.0], [6.0, 2.0, -1.0], [5.0, 3.0, -1.0]];
    let mapping = PointMapping::new(
        (0..3).map(|index| PointMatch {
            reference: index,
            model: index,
        }),
        reference.len(),
        model.len(),
    )
    .unwrap_or_else(|error| panic!("mapping failed: {error}"));
    let alignment = align_mapping(&reference, &model, &mapping)
        .unwrap_or_else(|error| panic!("alignment failed: {error}"));
    let measurement = measure_mapping(&reference, &model, &mapping, alignment);
    let verdict = decide_rmsd(&measurement, 1.0e-5);
    assert_eq!(measurement.distances.len(), 3);
    assert!(measurement.rmsd < 1.0e-6);
    assert!(verdict.passed);
}

#[test]
fn mapping_refuses_reusing_one_model_atom() {
    let result = PointMapping::new(
        [
            PointMatch {
                reference: 0,
                model: 0,
            },
            PointMatch {
                reference: 1,
                model: 0,
            },
        ],
        2,
        2,
    );
    assert!(matches!(result, Err(CompareError::InvalidMapping)));
}
