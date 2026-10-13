use super::*;

#[test]
fn the_default_source_is_bondi_and_says_so() {
    let source = RadiusSource::default();
    assert_eq!(source.radius(Element::CARBON), Some(1.70));
    assert_eq!(source.identity(), RadiusSet::Bondi.identity());
    assert_eq!(source.set(), Some(RadiusSet::Bondi));
}

#[test]
fn a_custom_source_reports_the_callers_identity() {
    let Ok(table) = CustomRadii::parse("name: mine\nversion: 9\nsource: me\nlicence: MIT\nC 2.0\n")
    else {
        panic!("table parses");
    };
    let source = RadiusSource::from(table);
    assert_eq!(source.radius(Element::CARBON), Some(2.0));
    assert_eq!(source.identity().name, "mine");
    assert_eq!(source.set(), None);
}
