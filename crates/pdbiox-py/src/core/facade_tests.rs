use super::{PyBondInference, PyProteinAlphaTrace, PySideChainTorsionReport};

#[test]
fn facade_defaults_are_exposed_without_reinterpreting_native_values() {
    let options = PyBondInference::standard();
    let native = pdbiox::BondInference::default();
    assert!((options.scale() - native.scale).abs() < f32::EPSILON);
    assert!((options.lower_bound() - native.lower_bound).abs() < f32::EPSILON);
    assert_eq!(
        options.exclude_across_chains(),
        native.exclude_across_chains
    );
    assert_eq!(options.respect_existing(), native.respect_existing);
}

#[test]
fn projection_types_keep_their_python_owned_shapes() {
    let trace = PyProteinAlphaTrace {
        chain: 2,
        positions: vec![Some([1.0, 2.0, 3.0]), None],
    };
    assert_eq!(trace.chain, 2);
    assert_eq!(trace.positions.len(), 2);

    let report = PySideChainTorsionReport {
        records: Vec::new(),
        findings: vec!["finding".to_owned()],
        dictionary_version: "ccd".to_owned(),
    };
    assert_eq!(report.findings, ["finding"]);
    assert_eq!(report.dictionary_version, "ccd");
}
