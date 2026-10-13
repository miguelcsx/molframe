use super::*;

const HEADER: &str = "name: lab\nversion: 2\nsource: bench\nlicence: CC0-1.0\n";

fn default_role(residue: &str, atom: &str) -> Option<u8> {
    HeavyRoleTable::builtin().role(residue, atom, false, true)
}

#[test]
fn default_text_parses_and_names_its_source() {
    let parsed = build(HEAVY_ROLES_DEFAULT_TEXT, Some(HEAVY_ROLES_DEFAULT)).expect("parses");
    assert_eq!(parsed.identity().name, "hbplus-baker-hubbard-heavy-roles");
    assert!(!parsed.exact.is_empty());
    assert_eq!(HeavyRoleTable::builtin(), &parsed);
}

#[test]
fn default_rows_reproduce_the_documented_conventions() {
    assert_eq!(default_role("ALA", "N"), Some(DONOR));
    assert_eq!(default_role("ALA", "OXT"), Some(ACCEPTOR));
    assert_eq!(default_role("ALA", "CB"), Some(0));
    assert_eq!(default_role("SER", "OG"), Some(BOTH));
    assert_eq!(default_role("HID", "ND1"), Some(DONOR));
    assert_eq!(default_role("HID", "NE2"), Some(ACCEPTOR));
    assert_eq!(default_role("HID", "N"), Some(DONOR));
    assert_eq!(default_role("HOH", "O"), Some(BOTH));
    assert_eq!(default_role("HOH", "H1"), Some(0));
    assert_eq!(default_role("DG", "N1"), Some(DONOR));
    assert_eq!(default_role("A", "O2'"), Some(BOTH));
    assert_eq!(default_role("ZZZ", "N"), None);
}

#[test]
fn weak_and_chain_start_conditions() {
    let table = HeavyRoleTable::builtin();
    assert_eq!(table.role("PRO", "N", false, false), Some(0));
    assert_eq!(table.role("PRO", "N", true, false), Some(DONOR));
    assert_eq!(table.role("CYS", "SG", false, false), Some(0));
    assert_eq!(table.role("CYS", "SG", false, true), Some(BOTH));
    assert_eq!(table.role("CYX", "SG", false, false), Some(0));
    assert_eq!(table.role("CYX", "SG", false, true), Some(ACCEPTOR));
    assert_eq!(table.role("MET", "SD", false, true), Some(ACCEPTOR));
}

#[test]
fn a_supplied_table_replaces_the_defaults() {
    let text = format!("{HEADER}amino: ALA\nwater: HOH\n@amino N acceptor\n@water O* donor\n");
    let table = HeavyRoleTable::parse(&text).expect("parses");
    assert_eq!(table.identity().name, "lab");
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
        let error = HeavyRoleTable::parse(text).expect_err(needle).to_string();
        assert!(error.contains(needle), "{error} lacks {needle}");
    }
    let error = HeavyRoleTable::parse(&format!("{HEADER}ALA N donor\nALA O maybe\n"))
        .expect_err("bad role");
    assert_eq!(error.line(), 6);
}
