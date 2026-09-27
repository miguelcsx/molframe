use crate::predicate_pattern::{NumericPattern, parse_residue, split_range};
use molframe_core::io::{InputBuffer, ReadOptions};

#[test]
fn signed_numbers_and_insertion_codes_split_without_losing_the_sign() {
    let negative = match parse_residue("-12A") {
        Ok(value) => value,
        Err(error) => panic!("parse failed: {error}"),
    };
    assert_eq!(negative.number, -12);
    assert_eq!(negative.insertion.as_ref(), "A");
    assert_eq!(split_range("-12A--10B"), Some(("-12A", "-10B")));
}

#[test]
fn numeric_ranges_are_inclusive_and_order_independent() {
    let pattern = match NumericPattern::parse("5:2") {
        Ok(pattern) => pattern,
        Err(error) => panic!("parse failed: {error}"),
    };
    assert!(pattern.matches(2.0));
    assert!(pattern.matches(5.0));
    assert!(!pattern.matches(5.1));
}

/// Author and label names as a deposited mmCIF carries them, with an iron
/// typed in upper case.
const DEPOSITED: &str = "data_d\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n\
_atom_site.auth_atom_id\n\
ATOM 1 N N VAL A 1 0 0 0 1 VAL A N\n\
ATOM 2 C CA VAL A 1 1 0 0 1 VAL A CA\n\
HETATM 3 FE FE HEM B . 3 0 0 142 HEM A FE\n";

fn deposited_count(source: &str) -> u64 {
    let input = InputBuffer::from_bytes(DEPOSITED.as_bytes().to_vec());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let query = match crate::Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("{source} compiles: {findings:?}"),
    };
    match query.evaluate(
        &structure,
        &molframe_core::contract::AnalysisPolicy::default(),
        &crate::Groups::new(),
        None,
    ) {
        Ok(evaluation) => evaluation.selection.len(),
        Err(findings) => panic!("{source} evaluates: {findings:?}"),
    }
}

#[test]
fn atom_names_match_under_the_default_identifier_policy() {
    assert_eq!(deposited_count("name CA"), 1);
    assert_eq!(deposited_count("name N CA"), 2);
    assert_eq!(deposited_count("name C*"), 1);
    assert_eq!(deposited_count("label_name CA"), 1);
    assert_eq!(deposited_count("auth_name CA"), 1);
}

#[test]
fn an_element_matches_whatever_case_the_query_and_file_use() {
    for query in ["element Fe", "element FE", "element fe"] {
        assert_eq!(deposited_count(query), 1, "{query}");
    }
    assert_eq!(deposited_count("element C"), 1);
    assert_eq!(deposited_count("element N C"), 2);
}
