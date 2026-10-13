use super::*;
use crate::builtin_validation::fixtures::reference;
use molframe_chem::TableError;

#[test]
fn a_given_table_resolves_with_its_identity() -> Result<(), BuiltinError> {
    let set = reference();
    assert_eq!(set.geometry()?.identity.name, "synthetic-geometry");
    assert_eq!(set.rama()?.identity.name, "synthetic-rama");
    assert_eq!(set.chi()?.table.len(), 1);
    Ok(())
}

#[test]
fn a_set_with_nothing_selected_names_how_to_supply_each_table() {
    let directory = std::env::temp_dir().join("molframe-empty-reference-tests");
    let set = ReferenceSet::with_store(Arc::new(TableStore::new(directory)));
    let missing = |result: Result<(), BuiltinError>| match result {
        Err(BuiltinError::Table(TableError::Missing { option, .. })) => Some(option),
        _ => None,
    };
    assert_eq!(missing(set.geometry().map(drop)), Some("--geometry-table"));
    assert_eq!(missing(set.rama().map(drop)), Some("--rama-table"));
    assert_eq!(missing(set.rotamers().map(drop)), Some("--rotamer-table"));
    assert_eq!(missing(set.chi().map(drop)), Some("--chi-table"));
    assert_eq!(missing(set.nucleic().map(drop)), Some("--nucleic-table"));
}
