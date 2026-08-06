use super::*;

#[test]
fn success_is_zero_and_the_codes_are_the_documented_ones() {
    assert_eq!(Exit::Success.code(), 0);
    assert_eq!(Exit::Usage.code(), 2);
    assert_eq!(Exit::Parse.code(), 4);
    assert_eq!(Exit::Refused.code(), 7);
    assert_eq!(Exit::Indeterminate.code(), 9);
    assert_eq!(Exit::Resource.code(), 10);
}

#[test]
fn no_findings_means_success() {
    assert_eq!(Exit::of(&[]), Exit::Success);
}

#[test]
fn a_refused_conversion_outranks_everything_else() {
    let findings = [
        Diagnostic::new(Code::W2001),
        Diagnostic::new(Code::E4102),
        Diagnostic::new(Code::E1103),
    ];
    assert_eq!(Exit::of(&findings), Exit::Refused);
}

#[test]
fn each_class_of_finding_maps_onto_its_own_code() {
    assert_eq!(Exit::of(&[Diagnostic::new(Code::E1103)]), Exit::Parse);
    assert_eq!(Exit::of(&[Diagnostic::new(Code::E2001)]), Exit::Schema);
    assert_eq!(Exit::of(&[Diagnostic::new(Code::E3005)]), Exit::Consistency);
    assert_eq!(Exit::of(&[Diagnostic::new(Code::E6002)]), Exit::Policy);
    assert_eq!(Exit::of(&[Diagnostic::new(Code::E1901)]), Exit::Resource);
    assert_eq!(Exit::of(&[Diagnostic::new(Code::E1001)]), Exit::Input);
}

#[test]
fn the_worst_of_several_findings_decides() {
    let findings = [Diagnostic::new(Code::W2001), Diagnostic::new(Code::E3005)];
    assert_eq!(Exit::of(&findings), Exit::Consistency);
}
