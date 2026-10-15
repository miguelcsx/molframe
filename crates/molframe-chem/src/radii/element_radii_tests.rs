use super::*;

const HEADER: &str = "name: lab\nversion: 2\nsource: survey\nlicence: CC0-1.0\n";

fn table(body: &str) -> Result<Provided<VdwRadii>, TableTextError> {
    VdwRadii::decode(
        format!("{HEADER}{body}").as_bytes(),
        TableIdentity::supplied("x", "x", "x", "x"),
    )
}

#[test]
fn a_table_declares_its_identity_and_its_radii() {
    let Ok(provided) = table("# note\nC 1.75\nn 1.6\n") else {
        panic!("a well-formed table parses");
    };
    assert_eq!(provided.identity.name, "lab");
    assert_eq!(provided.identity.licence, "CC0-1.0");
    assert_eq!(provided.table.radius(Element::CARBON), Some(1.75));
    assert_eq!(provided.table.radius(Element::NITROGEN), Some(1.6));
    assert_eq!(provided.table.radius(Element::OXYGEN), None);
}

#[test]
fn the_two_roles_share_one_reader_and_differ_in_their_names() {
    assert_eq!(VdwRadii::KIND, "vdw-radii");
    assert_eq!(CovalentRadii::KIND, "covalent-radii");
    assert_ne!(VdwRadii::EXTENSION, CovalentRadii::EXTENSION);
    let Ok(provided) = CovalentRadii::decode(
        format!("{HEADER}C 0.76\n").as_bytes(),
        TableIdentity::supplied("x", "x", "x", "x"),
    ) else {
        panic!("a covalent table parses");
    };
    assert_eq!(provided.table.radius(Element::CARBON), Some(0.76));
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
    let Err(error) = VdwRadii::decode(
        b"name: a\nversion: 1\nsource: s\nC 1.7\n",
        TableIdentity::supplied("x", "x", "x", "x"),
    ) else {
        panic!("a table without a licence is refused");
    };
    assert!(error.to_string().contains("licence"));
    assert_eq!(error.line(), 0);
    assert!(table("").is_err());
}
