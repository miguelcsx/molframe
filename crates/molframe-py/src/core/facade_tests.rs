use super::{PyBondInference, PyProteinAlphaTrace, PySideChainTorsionReport};

#[test]
fn facade_defaults_are_exposed_without_reinterpreting_native_values() {
    let options = PyBondInference::standard();
    let native = molframe::BondInference::default();
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
        findings: vec![
            molframe::Diagnostic::new(molframe::Code::W2001)
                .at_row(7)
                .into(),
        ],
        dictionary_version: "ccd".to_owned(),
    };
    // The row is a structured detail the text form used to flatten away.
    assert_eq!(report.findings[0].inner.row(), Some(7));
    assert_eq!(report.dictionary_version, "ccd");
}
