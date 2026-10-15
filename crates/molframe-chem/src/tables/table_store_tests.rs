use super::*;
use crate::{Provided, TableIdentity, TableTextError};

#[derive(Debug, PartialEq)]
struct Words(Vec<String>);

impl TableKind for Words {
    const KIND: &'static str = "words";
    const EXTENSION: &'static str = "txt";

    fn decode(bytes: &[u8], fallback: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let text = crate::text_of("words table", bytes)?;
        Ok(Provided {
            identity: fallback,
            table: Self(text.split_whitespace().map(str::to_owned).collect()),
        })
    }
}

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
