use super::*;

#[test]
fn codes_render_with_a_stable_prefix_and_four_digits() {
    assert_eq!(Code::E1001.to_string(), "MOLFRAME-E1001");
    assert_eq!(Code::W3011.to_string(), "MOLFRAME-W3011");
}

#[test]
fn strictness_decides_which_severities_are_errors() {
    assert!(!Severity::Loose.is_error(Strictness::Strict));
    assert!(Severity::Strict.is_error(Strictness::Strict));
    assert!(!Severity::Strict.is_error(Strictness::Medium));
    assert!(Severity::Invalidating.is_error(Strictness::Medium));
    assert!(!Severity::Invalidating.is_error(Strictness::Loose));
    assert!(Severity::Breaking.is_error(Strictness::Loose));
}

#[test]
fn the_leading_digit_names_the_class() {
    assert_eq!(Code::E1103.class(), Class::Syntax);
    assert_eq!(Code::E2001.class(), Class::Schema);
    assert_eq!(Code::W3011.class(), Class::Consistency);
    assert_eq!(Code::E4102.class(), Class::Conversion);
    assert_eq!(Code::E6001.class(), Class::Policy);
    assert_eq!(Code::E9001.class(), Class::Internal);
}

#[test]
fn an_unregistered_code_reads_as_breaking_rather_than_as_harmless() {
    let invented = Code::new(Kind::Error, 1);
    assert!(!invented.is_registered());
    assert_eq!(invented.severity(), Severity::Breaking);
}
