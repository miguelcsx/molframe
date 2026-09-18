use super::*;
use crate::contract::{AnalysisPolicy, SourceRef};

#[test]
fn json_names_every_policy_field_in_stable_order_and_escapes_the_source() {
    let record = Provenance::new(&AnalysisPolicy::default())
        .with_source(SourceRef::Url("https://example.test/a\"b".into()));
    let first = record.to_json();
    let second = record.to_json();
    assert_eq!(first, second);
    assert!(first.contains("\"input_source\":\"https://example.test/a\\\"b\""));
    assert!(first.contains("\"policy.atom_equivalence\":\"Ccd\""));
    assert!(first.contains("\"policy.float_tolerance.absolute\":"));
}

#[test]
fn mmcif_is_a_standalone_ordered_category_and_quotes_spaces() {
    let record = Provenance::new(&AnalysisPolicy::default())
        .with_source(SourceRef::path("a path/input.cif"));
    let cif = record.to_mmcif();
    assert!(cif.starts_with("data_molframe_provenance\n#\nloop_\n"));
    assert!(cif.contains("input_source 'a path/input.cif'"));
    assert!(cif.contains("policy.missing_atoms Report"));
    assert!(cif.ends_with("#\n"));
}
