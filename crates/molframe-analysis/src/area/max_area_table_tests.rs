use super::MaxAreaTable;
use molframe_chem::{Provided, TableIdentity, TableKind};

const HEADER: &str = "name: lab\nversion: 1\nsource: our ensemble\nlicence: CC0\n";

fn decode(text: &str) -> Result<Provided<MaxAreaTable>, molframe_chem::TableTextError> {
    MaxAreaTable::decode(
        text.as_bytes(),
        TableIdentity::supplied("fallback", "0", "none", "none"),
    )
}

#[test]
fn a_caller_table_declares_its_identity_and_values() {
    let Ok(provided) = decode(&format!("{HEADER}# c\nala 130.5\nGLY 105\n")) else {
        panic!("a valid table");
    };
    assert!(provided.identity.to_string().contains("lab 1"));
    assert_eq!(provided.table.max_area("ALA"), Some(130.5));
    assert_eq!(provided.table.max_area("TRP"), None);
    assert_eq!(provided.table.relative("GLY", 52.5), Some(0.5));
}

#[test]
fn a_variant_takes_its_parents_value() {
    let Ok(provided) = decode(&format!("{HEADER}HIS 200\n")) else {
        panic!("a valid table");
    };
    assert_eq!(provided.table.max_area("HIE"), Some(200.0));
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
        let Err(error) = decode(text) else {
            panic!("{text} must be refused");
        };
        assert_eq!(error.line(), line, "{text}");
        assert!(error.to_string().contains(fragment), "{error}");
    }
}

/// The shipped Tien table, when the repository's data directory is present.
#[test]
fn the_shipped_table_reproduces_the_published_maxima() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/tien-2013-max-area.maxarea"
    );
    let Ok(bytes) = std::fs::read(path) else {
        return;
    };
    let Ok(provided) = MaxAreaTable::decode(&bytes, TableIdentity::supplied("", "", "", "")) else {
        panic!("the shipped table must decode");
    };
    assert_eq!(provided.identity.name, "tien-2013-theoretical");
    assert_eq!(provided.table.max_area("GLY"), Some(104.0));
    assert_eq!(provided.table.max_area("HID"), Some(224.0));
    assert_eq!(provided.table.max_area("XYZ"), None);
}
