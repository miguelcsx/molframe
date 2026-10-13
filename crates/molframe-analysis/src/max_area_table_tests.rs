use super::MaxAreaTable;

const HEADER: &str = "name: lab\nversion: 1\nsource: our ensemble\nlicence: CC0\n";

#[test]
fn the_default_is_the_tien_table_and_names_itself() {
    let table = MaxAreaTable::default();
    assert_eq!(table.identity().name, "tien-2013-theoretical");
    assert_eq!(table.max_area("GLY"), Some(104.0));
    assert_eq!(table.max_area("HID"), Some(224.0));
    assert_eq!(table.max_area("XYZ"), None);
    assert_eq!(table.relative("GLY", 52.0), Some(0.5));
}

#[test]
fn a_caller_table_replaces_the_values_and_declares_its_identity() {
    let table = MaxAreaTable::parse(&format!("{HEADER}# c\nala 130.5\nGLY 105\n")).expect("valid");
    assert!(table.identity().to_string().contains("lab 1"));
    assert_eq!(table.max_area("ALA"), Some(130.5));
    assert_eq!(table.max_area("TRP"), None);
}

#[test]
fn a_variant_takes_its_parents_value_from_the_callers_table() {
    let table = MaxAreaTable::parse(&format!("{HEADER}HIS 200\n")).expect("valid");
    assert_eq!(table.max_area("HIE"), Some(200.0));
}

#[test]
fn every_problem_is_refused_with_its_line() {
    let cases = [
        ("name: a\nversion: 1\nsource: s\nALA 1\n", 0, "licence"),
        (&format!("{HEADER}ALA -1\n"), 5, "positive"),
        (&format!("{HEADER}ALA x\n"), 5, "not a number"),
        (&format!("{HEADER}ALA 1\nALA 2\n"), 6, "twice"),
        (&format!("{HEADER}ALA 1 2\n"), 5, "two words"),
        (
            &format!("{HEADER}colour: red\nALA 1\n"),
            5,
            "unknown header",
        ),
        (HEADER, 0, "no residue"),
        (&format!("{HEADER}ALA 9999\n"), 5, "at most"),
    ];
    for (text, line, fragment) in cases {
        let error = MaxAreaTable::parse(text).expect_err(text);
        assert_eq!(error.line(), line, "{text}");
        assert!(error.to_string().contains(fragment), "{error}");
    }
}
