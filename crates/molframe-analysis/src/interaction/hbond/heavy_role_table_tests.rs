use super::*;
use molframe_chem::{Provided, TableIdentity};

const HEADER: &str = "name: lab\nversion: 2\nsource: bench\nlicence: CC0-1.0\n";

fn decode(text: &str) -> Result<Provided<HeavyRoleTable>, TableTextError> {
    HeavyRoleTable::decode(
        text.as_bytes(),
        TableIdentity::supplied("fallback", "0", "none", "none"),
    )
}

const SYNTHETIC: &str = "amino: ALA PRO CYS HIS SER\nnucleic: A\nwater: HOH\n\
HID ND1 donor\nHID NE2 acceptor\nCYX SG acceptor weak\n@water O* both\n\
PRO N donor chain-start\nPRO N none\n@amino N donor\n@amino O,OXT acceptor\n\
SER OG both\nCYS SG both weak\n@nucleic O2' both\n";

fn synthetic() -> HeavyRoleTable {
    match decode(&format!("{HEADER}{SYNTHETIC}")) {
        Ok(provided) => provided.table,
        Err(error) => panic!("{error}"),
    }
}

#[test]
fn rows_reproduce_the_conventions_they_declare() {
    let table = synthetic();
    let role = |residue: &str, atom: &str| table.role(residue, atom, false, true);
    assert_eq!(role("ALA", "N"), Some(DONOR));
    assert_eq!(role("ALA", "OXT"), Some(ACCEPTOR));
    assert_eq!(role("ALA", "CB"), Some(0));
    assert_eq!(role("SER", "OG"), Some(BOTH));
    assert_eq!(role("HID", "ND1"), Some(DONOR));
    assert_eq!(role("HID", "NE2"), Some(ACCEPTOR));
    assert_eq!(role("HOH", "O"), Some(BOTH));
    assert_eq!(role("HOH", "H1"), Some(0));
    assert_eq!(role("A", "O2'"), Some(BOTH));
    assert_eq!(role("ZZZ", "N"), None);
}

#[test]
fn weak_and_chain_start_conditions() {
    let table = synthetic();
    assert_eq!(table.role("PRO", "N", false, false), Some(0));
    assert_eq!(table.role("PRO", "N", true, false), Some(DONOR));
    assert_eq!(table.role("CYS", "SG", false, false), Some(0));
    assert_eq!(table.role("CYS", "SG", false, true), Some(BOTH));
    assert_eq!(table.role("CYX", "SG", false, false), Some(0));
    assert_eq!(table.role("CYX", "SG", false, true), Some(ACCEPTOR));
}

/// The shipped table decodes and names its source, when the data directory is present.
#[test]
fn the_shipped_table_decodes() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/heavy-atom-roles.roles"
    );
    let Ok(bytes) = std::fs::read(path) else {
        return;
    };
    let Ok(provided) = HeavyRoleTable::decode(&bytes, TableIdentity::supplied("", "", "", ""))
    else {
        panic!("the shipped table must decode");
    };
    assert_eq!(provided.identity.name, "hbplus-baker-hubbard-heavy-roles");
    assert_eq!(provided.table.role("HID", "N", false, true), Some(DONOR));
}

#[test]
fn a_supplied_table_replaces_the_defaults() {
    let text = format!("{HEADER}amino: ALA\nwater: HOH\n@amino N acceptor\n@water O* donor\n");
    let Ok(provided) = decode(&text) else {
        panic!("a valid table");
    };
    let table = provided.table;
    assert_eq!(provided.identity.name, "lab");
    assert_eq!(table.role("ALA", "N", false, false), Some(ACCEPTOR));
    assert_eq!(table.role("HOH", "OW", false, false), Some(DONOR));
    assert_eq!(table.role("SER", "OG", false, false), None);
}

#[test]
fn malformed_tables_name_the_line() {
    let cases = [
        ("amino: ALA\nA N donor\n", "name"),
        (
            &format!("{HEADER}amino: ALA\nALA N sometimes\n"),
            "unknown role",
        ),
        (&format!("{HEADER}ALA N donor extra tail\n"), "condition"),
        (&format!("{HEADER}@amino N donor\n"), "before its class"),
        (
            &format!("{HEADER}amino: ALA\namino: GLY\nALA N donor\n"),
            "twice",
        ),
        (
            &format!("{HEADER}colour: red\nALA N donor\n"),
            "unknown header",
        ),
        (&format!("{HEADER}ALA N\n"), "expected"),
        (HEADER, "no row"),
    ];
    for (text, needle) in cases {
        let error = decode(text).expect_err(needle).to_string();
        assert!(error.contains(needle), "{error} lacks {needle}");
    }
    let error = decode(&format!("{HEADER}ALA N donor\nALA O maybe\n")).expect_err("bad role");
    assert_eq!(error.line(), 6);
}
