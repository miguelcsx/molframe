use super::*;
use molframe::{Code, Diagnostic};

#[test]
fn every_registered_code_names_a_class_this_module_can_raise() {
    for code in Code::registered() {
        assert!(
            CLASS_NAMES.contains(&class_name(code)),
            "{code} maps to {}, which is not a declared class",
            class_name(code)
        );
    }
}

#[test]
fn the_exact_codes_win_over_their_class() {
    assert_eq!(class_name(Code::E7001), "MemoryBudgetError");
    assert_eq!(class_name(Code::E1902), "MemoryBudgetError");
    assert_eq!(class_name(Code::E1904), "Cancelled");
    assert_eq!(class_name(Code::E1901), "ResourceError");
    assert_eq!(class_name(Code::E7101), "MolframeIOError");
    assert_eq!(class_name(Code::E7901), "MolframeIOError");
    assert_eq!(class_name(Code::E1103), "ParseError");
    assert_eq!(class_name(Code::E5101), "GeometryError");
    assert_eq!(class_name(Code::E6101), "PolicyError");
    assert_eq!(class_name(Code::E9001), "InternalError");
}

#[test]
fn a_diagnostic_becomes_an_exception_with_its_structured_fields() {
    Python::initialize();
    Python::attach(|py| {
        install_package_for_tests(py);
        let diagnostic = Diagnostic::new(Code::E5102).with_message("lengths differ");
        let error = from_diagnostic(&diagnostic);
        let value = error.value(py);
        let classes = errors_module(py).expect("molframe.errors imports");
        let geometry = classes.getattr("GeometryError").expect("class exists");
        assert!(value.is_instance(&geometry).expect("isinstance works"));
        let code: String = value.getattr("code").expect("code").extract().expect("str");
        assert_eq!(code, "MOLFRAME-E5102");
        let message: String = value
            .getattr("message")
            .expect("message")
            .extract()
            .expect("str");
        assert_eq!(message, "lengths differ");
        let remedy: String = value
            .getattr("remedy")
            .expect("remedy")
            .extract()
            .expect("str");
        assert!(!remedy.is_empty());
        let findings = value.getattr("findings").expect("findings");
        assert_eq!(findings.len().expect("a tuple"), 1);
    });
}

#[test]
fn a_memory_budget_failure_is_also_a_memory_error() {
    Python::initialize();
    Python::attach(|py| {
        install_package_for_tests(py);
        let error = from_diagnostic(&Diagnostic::new(Code::E7001));
        assert!(error.is_instance_of::<pyo3::exceptions::PyMemoryError>(py));
    });
}

#[test]
fn a_bad_value_is_also_a_value_error_and_a_missing_key_a_key_error() {
    Python::initialize();
    Python::attach(|py| {
        install_package_for_tests(py);
        assert!(value("no").is_instance_of::<pyo3::exceptions::PyValueError>(py));
        assert!(key("no").is_instance_of::<pyo3::exceptions::PyKeyError>(py));
        assert!(index("no").is_instance_of::<pyo3::exceptions::PyIndexError>(py));
        assert!(type_error("no").is_instance_of::<pyo3::exceptions::PyTypeError>(py));
    });
}
