use super::*;
use crate::TableKind;

fn given(name: &str, carbon: f32) -> Provided<VdwRadii> {
    let text = format!("name: {name}\nversion: 9\nsource: me\nlicence: MIT\nC {carbon}\n");
    let Ok(provided) =
        VdwRadii::decode(text.as_bytes(), TableIdentity::supplied("x", "x", "x", "x"))
    else {
        panic!("table parses");
    };
    provided
}

#[test]
fn a_source_answers_from_its_table_and_reports_its_identity() {
    let source = RadiusSource::from(given("mine", 2.0));
    assert_eq!(source.radius(Element::CARBON), Some(2.0));
    assert_eq!(source.radius(Element::OXYGEN), None);
    assert_eq!(source.identity().name, "mine");
}

#[test]
fn a_request_resolves_through_the_store_and_never_picks_silently() {
    // No compiled-in layer either: nothing may answer for the profile.
    let store = TableStore::new(std::env::temp_dir().join("molframe-radius-source-empty"))
        .with_bundled(None);
    let request = TableRequest::given(given("held", 1.5));
    let Ok(source) = RadiusSource::resolve(&request, &store) else {
        panic!("a given table resolves");
    };
    assert_eq!(source.identity().name, "held");
    let Err(error) = RadiusSource::resolve(&TableRequest::profile(), &store) else {
        panic!("no profile and no request is an error");
    };
    assert!(error.to_string().contains("--radii"), "{error}");
}
