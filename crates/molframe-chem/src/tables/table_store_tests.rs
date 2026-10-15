use super::*;
use crate::test_data::Words;

fn store(name: &str) -> TableStore {
    let dir = std::env::temp_dir().join(format!("molframe-store-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::create_dir_all(&dir);
    TableStore::new(dir)
}

fn put(store: &TableStore, file: &str, text: &str) {
    let _ = std::fs::write(store.directory().join(file), text);
}

#[test]
fn nothing_requested_and_no_profile_is_a_typed_missing() {
    let store = store("none");
    let Err(error) = store.select::<Words>(None, "--words") else {
        panic!("nothing selects a table");
    };
    assert!(matches!(error, TableError::Missing { .. }));
    assert!(error.to_string().contains("--words"));
}

#[test]
fn the_profile_names_the_entry_a_kind_uses() {
    let store = store("profile");
    put(&store, "lab.txt", "a b");
    put(&store, "profile.txt", "# comment\nother: x\nwords: lab\n");
    let Ok(chosen) = store.select::<Words>(None, "--words") else {
        panic!("the profile selects lab");
    };
    assert_eq!(chosen.source, Source::Cache);
    assert_eq!(chosen.provided.table, Words(vec!["a".into(), "b".into()]));
}

#[test]
fn a_profile_without_the_kind_is_missing() {
    let store = store("other");
    put(&store, "profile.txt", "geometry: g\n");
    assert!(matches!(
        store.select::<Words>(None, "--words"),
        Err(TableError::Missing { .. })
    ));
}

#[test]
fn a_named_entry_the_cache_lacks_is_missing_not_unreadable() {
    let store = store("lacks");
    assert!(matches!(
        store.select::<Words>(Some("ghost"), "--words"),
        Err(TableError::Missing { .. })
    ));
}

#[test]
fn a_path_request_is_a_user_file_and_beats_the_profile() {
    let store = store("path");
    put(&store, "lab.txt", "cache");
    put(&store, "profile.txt", "words: lab\n");
    let file = store.directory().join("mine.dat");
    let _ = std::fs::write(&file, "user");
    let Ok(chosen) = store.select::<Words>(file.to_str(), "--words") else {
        panic!("the path is read");
    };
    assert_eq!(chosen.source, Source::User);
    assert_eq!(chosen.provided.table, Words(vec!["user".into()]));
}

#[test]
fn a_table_is_decoded_once_and_read_again_after_an_edit() {
    let store = store("memo");
    put(&store, "lab.txt", "a");
    let (Ok(first), Ok(second)) = (
        store.select::<Words>(Some("lab"), "--words"),
        store.select::<Words>(Some("lab"), "--words"),
    ) else {
        panic!("lab loads");
    };
    assert!(Arc::ptr_eq(&first.provided, &second.provided));
    put(&store, "lab.txt", "a b c");
    let Ok(third) = store.select::<Words>(Some("lab"), "--words") else {
        panic!("lab reloads");
    };
    assert_eq!(third.provided.table.0.len(), 3);
}

/// A defaults layer standing in for the compiled-in one, so precedence is
/// tested without the feature.
fn defaults(profile: &str, entries: &[(&str, &'static str)]) -> Arc<MemoryTables> {
    let mut tables = MemoryTables::new("defaults");
    tables.set_profile(profile);
    for (name, text) in entries {
        tables.insert(name, "txt", TableBytes::from(*text));
    }
    Arc::new(tables)
}

fn words(chosen: &Selected<Words>) -> String {
    chosen.provided.table.0.join(" ")
}

#[test]
fn the_defaults_answer_when_there_is_no_profile_or_cache() {
    let store =
        store("defaults").with_bundled(Some(defaults("words: lab\n", &[("lab", "fallback")])));
    let Ok(Some(chosen)) = store.select_optional::<Words>("--words") else {
        panic!("the defaults select lab");
    };
    assert_eq!(chosen.source, Source::Memory);
    assert_eq!(words(&chosen), "fallback");
}

#[test]
fn a_profile_on_disk_wins_over_the_defaults() {
    let store = store("disk-wins").with_bundled(Some(defaults(
        "words: lab\n",
        &[("lab", "fallback"), ("mine", "unused")],
    )));
    put(&store, "profile.txt", "words: mine\n");
    put(&store, "mine.txt", "disk");
    let Ok(chosen) = store.select::<Words>(None, "--words") else {
        panic!("the disk profile selects mine");
    };
    assert_eq!(chosen.source, Source::Cache);
    assert_eq!(words(&chosen), "disk");
}

#[test]
fn a_cache_file_wins_over_a_default_of_the_same_name() {
    let store =
        store("disk-entry").with_bundled(Some(defaults("words: lab\n", &[("lab", "fallback")])));
    put(&store, "lab.txt", "disk");
    let Ok(chosen) = store.select::<Words>(None, "--words") else {
        panic!("lab is selected");
    };
    assert_eq!(chosen.source, Source::Cache);
    assert_eq!(words(&chosen), "disk");
}

#[test]
fn a_disk_profile_without_the_kind_defers_to_the_defaults() {
    let store = store("defer").with_bundled(Some(defaults("words: lab\n", &[("lab", "fallback")])));
    put(&store, "profile.txt", "geometry: g\n");
    let Ok(chosen) = store.select::<Words>(None, "--words") else {
        panic!("the defaults still name lab");
    };
    assert_eq!(words(&chosen), "fallback");
}

#[test]
fn a_memory_registration_wins_over_the_disk_and_the_defaults() {
    let store =
        store("memory-wins").with_bundled(Some(defaults("words: lab\n", &[("lab", "fallback")])));
    put(&store, "profile.txt", "words: lab\n");
    put(&store, "lab.txt", "disk");
    store.insert_table::<Words>("lab", "memory");
    let Ok(chosen) = store.select::<Words>(None, "--words") else {
        panic!("lab is selected");
    };
    assert_eq!(chosen.source, Source::Memory);
    assert_eq!(words(&chosen), "memory");
}

#[test]
fn a_memory_profile_wins_over_the_disk_profile() {
    let store = store("memory-profile");
    put(&store, "profile.txt", "words: lab\n");
    put(&store, "lab.txt", "disk");
    store.insert_profile("words: held\n");
    store.insert_entry("held", "txt", String::from("held"));
    let Ok(chosen) = store.select::<Words>(None, "--words") else {
        panic!("the memory profile selects held");
    };
    assert_eq!(words(&chosen), "held");
}

#[test]
fn two_selections_from_memory_share_one_decoded_table() {
    let store = store("shared");
    store.insert_profile("words: lab\n");
    store.insert_table::<Words>("lab", Arc::<[u8]>::from(b"a b".as_slice()));
    let (Ok(Some(first)), Ok(Some(second))) = (
        store.select_optional::<Words>("--words"),
        store.select_optional::<Words>("--words"),
    ) else {
        panic!("lab is selected");
    };
    assert!(Arc::ptr_eq(&first.provided, &second.provided));
}

#[test]
fn registering_an_entry_again_replaces_it() {
    let store = store("replace");
    store.insert_table::<Words>("lab", "a");
    store.insert_table::<Words>("lab", "a b c");
    let Ok(chosen) = store.select::<Words>(Some("lab"), "--words") else {
        panic!("lab is selected");
    };
    assert_eq!(chosen.provided.table.0.len(), 3);
}

#[test]
fn a_malformed_registration_is_reported_where_it_is_selected() {
    let store = store("malformed");
    store.insert_table::<Words>("lab", "bad");
    assert!(matches!(
        store.select::<Words>(Some("lab"), "--words"),
        Err(TableError::Malformed { .. })
    ));
}

#[test]
fn a_profile_that_cannot_be_read_is_unreadable_not_absent() {
    let store = store("unreadable");
    let profile = store.directory().join("profile-dir");
    let _ = std::fs::create_dir_all(&profile);
    let store = store.with_profile(profile);
    assert!(matches!(
        store.select_optional::<Words>("--words"),
        Err(TableError::Unreadable { .. })
    ));
}

#[test]
fn a_cache_directory_that_is_a_file_holds_nothing() {
    let store = store("not-a-dir");
    let file = store.directory().join("plain");
    let _ = std::fs::write(&file, "x");
    let store = TableStore::new(file.clone()).with_bundled(None);
    let store = store.with_profile(file.join("profile.txt"));
    assert!(matches!(
        store.select_optional::<Words>("--words"),
        Ok(None)
    ));
    assert!(matches!(
        store.select::<Words>(Some("lab"), "--words"),
        Err(TableError::Missing { .. })
    ));
}
