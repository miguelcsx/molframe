use super::*;
use crate::error::tests::install_package_for_tests;
use molframe::{Analysis, AnalysisPolicy, Coverage, Indeterminacy};

#[test]
fn an_indeterminate_analysis_has_no_value_and_asking_for_one_raises() {
    Python::initialize();
    Python::attach(|py| {
        install_package_for_tests(py);
        let policy = AnalysisPolicy::default();
        let refused = Analysis::<u32>::indeterminate(
            Indeterminacy::MissingInputs {
                missing: 3,
                ambiguous: 0,
            },
            Coverage::default(),
            &policy,
        );
        let envelope = PyAnalysis::new(&refused, None);
        assert_eq!(envelope.status(), "indeterminate");
        assert!(!envelope.is_determinate());
        assert_eq!(
            envelope.indeterminacy(),
            Some("the policy rejects 3 missing and 0 ambiguous inputs")
        );
        let Err(error) = envelope.value(py) else {
            panic!("an indeterminate analysis must not return a value");
        };
        let classes = py
            .import("molframe.errors")
            .expect("molframe.errors imports");
        let class = classes.getattr("IndeterminateError").expect("class exists");
        assert!(
            error
                .value(py)
                .is_instance(&class)
                .expect("isinstance works")
        );
    });
}

#[test]
fn a_determinate_analysis_returns_its_value() {
    Python::initialize();
    Python::attach(|py| {
        let answered = Analysis::complete(7_u32, Coverage::complete(1), &AnalysisPolicy::default());
        let value = 7_u32.into_pyobject(py).expect("an int").into_any().unbind();
        let envelope = PyAnalysis::new(&answered, Some(value));
        assert!(envelope.is_determinate());
        assert_eq!(envelope.indeterminacy(), None);
        let found: u32 = envelope
            .value(py)
            .expect("a value")
            .extract(py)
            .expect("an int");
        assert_eq!(found, 7);
    });
}
