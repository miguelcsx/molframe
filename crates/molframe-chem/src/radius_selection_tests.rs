use super::*;

fn store() -> TableStore {
    let data = std::path::PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data"));
    let profile = std::path::PathBuf::from(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../molframe-bench/data/profile.txt"
    ));
    TableStore::new(data).with_profile(profile)
}

fn name_of(selection: &RadiusSelection, set: RadiiSet) -> String {
    match selection.resolve(set, &store()) {
        Ok(source) => source.identity().name.to_string(),
        Err(error) => panic!("resolves: {error}"),
    }
}

#[test]
fn without_an_argument_the_policy_names_the_table_through_the_store() {
    let selection = RadiusSelection::choose(None, true);
    for set in [
        RadiiSet::Bondi,
        RadiiSet::Charmm,
        RadiiSet::Alvarez,
        RadiiSet::AmberUnited,
    ] {
        assert_eq!(name_of(&selection, set), set.name());
    }
}

#[test]
fn without_a_policy_the_argument_or_the_profile_decides() {
    assert_eq!(
        name_of(
            &RadiusSelection::choose(Some("alvarez"), false),
            RadiiSet::Charmm
        ),
        "alvarez"
    );
    assert_eq!(
        name_of(&RadiusSelection::choose(None, false), RadiiSet::Charmm),
        "bondi"
    );
}

#[test]
fn an_argument_equal_to_the_policy_is_accepted() {
    let selection = RadiusSelection::choose(Some("charmm"), true);
    assert_eq!(name_of(&selection, RadiiSet::Charmm), "charmm");
}

#[test]
fn an_argument_that_names_another_table_than_the_policy_is_a_conflict() {
    let selection = RadiusSelection::choose(Some("alvarez"), true);
    let Err(error) = selection.resolve(RadiiSet::Charmm, &store()) else {
        panic!("a conflict is an error");
    };
    assert!(matches!(
        error,
        TableError::Conflict { ref requested, ref policy, .. }
            if requested == "alvarez" && policy == "charmm"
    ));
    assert_eq!(
        molframe_core::Diagnostic::from(error).code(),
        molframe_core::Code::E7103
    );
}

#[test]
fn a_policy_table_the_store_lacks_is_missing_and_never_the_profile() {
    let empty = TableStore::new(std::env::temp_dir().join("molframe-no-such-radii-cache"));
    let Err(error) = RadiusSelection::Policy.resolve(RadiiSet::Charmm, &empty) else {
        panic!("a missing table is an error");
    };
    assert!(matches!(error, TableError::Missing { .. }));
}

#[test]
fn verify_ties_an_applied_table_to_the_policy_unless_it_stands_alone() {
    let applied = TableIdentity::supplied("bondi", "1964", "test", "CC0-1.0");
    assert!(
        RadiusSelection::Policy
            .verify(&applied, RadiiSet::Bondi)
            .is_ok()
    );
    assert!(
        RadiusSelection::Policy
            .verify(&applied, RadiiSet::Charmm)
            .is_err()
    );
    assert!(
        RadiusSelection::Table(TableRequest::Profile)
            .verify(&applied, RadiiSet::Charmm)
            .is_ok()
    );
}
