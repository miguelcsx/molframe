use super::*;

const HEADER: &str = "name: lab\nversion: 2\nsource: survey\nlicence: CC0-1.0\n";

fn table(body: &str) -> Result<CustomRadii, RadiusTableError> {
    CustomRadii::parse(&format!("{HEADER}{body}"))
}

#[test]
fn a_table_declares_its_identity_and_its_radii() {
    let Ok(radii) = table("# note\nC 1.75\nn 1.6\n") else {
        panic!("a well-formed table parses");
    };
    assert_eq!(radii.identity().name, "lab");
    assert_eq!(radii.identity().licence, "CC0-1.0");
    assert_eq!(radii.radius(Element::CARBON), Some(1.75));
    assert_eq!(radii.radius(Element::NITROGEN), Some(1.6));
    assert_eq!(radii.radius(Element::OXYGEN), None);
}

#[test]
fn a_base_fills_the_elements_the_table_leaves_out() {
    let Ok(radii) = table("base: bondi\nC 1.75\n") else {
        panic!("a layered table parses");
    };
    assert_eq!(radii.radius(Element::CARBON), Some(1.75));
    assert_eq!(radii.radius(Element::OXYGEN), Some(1.52));
    assert_eq!(radii.base(), Some(RadiusSet::Bondi));
}

#[test]
fn every_defect_is_refused_with_its_line() {
    let cases = [
        ("C 1.7 2.0\n", 5),
        ("Xx 1.7\n", 5),
        ("C abc\n", 5),
        ("C -1\n", 5),
        ("C 0\n", 5),
        ("C 99\n", 5),
        ("C 1.7\nC 1.8\n", 6),
        ("base: nothing\nC 1.7\n", 5),
        ("colour: red\nC 1.7\n", 5),
        ("name: again\nC 1.7\n", 5),
    ];
    for (body, line) in cases {
        let Err(error) = table(body) else {
            panic!("`{body}` must be refused");
        };
        assert_eq!(error.line(), line, "{body}: {error}");
    }
}

#[test]
fn missing_header_or_empty_table_is_refused_whole() {
    let Err(error) = CustomRadii::parse("name: a\nversion: 1\nsource: s\nC 1.7\n") else {
        panic!("a table without a licence is refused");
    };
    assert!(error.to_string().contains("licence"));
    assert_eq!(error.line(), 0);
    assert!(table("").is_err());
}
