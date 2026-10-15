use super::*;
use crate::test_data::Words;

fn directory(name: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("molframe-provider-{name}-{}", std::process::id()));
    let _ = std::fs::create_dir_all(&path);
    path
}

#[test]
fn a_user_file_is_decoded_once_and_borrowed_afterwards() {
    let dir = directory("user");
    let present = dir.join("t.txt");
    let _ = std::fs::write(&present, "x y");
    let file = UserFile::<Words>::new(present);
    let (Ok(first), Ok(second)) = (file.provide(), file.provide()) else {
        panic!("the file decodes");
    };
    assert!(core::ptr::eq(first, second));
    assert_eq!(first.table, Words(vec!["x".into(), "y".into()]));
    assert_eq!(first.identity.source, file_path_text(&dir.join("t.txt")));
    let Ok(shared) = file.shared() else {
        panic!("shared");
    };
    assert!(core::ptr::eq(Arc::as_ptr(&shared), first));
}

fn file_path_text(path: &std::path::Path) -> String {
    path.display().to_string()
}

#[test]
fn a_missing_user_file_is_an_error_not_a_fallthrough() {
    let dir = directory("missing");
    let absent = UserFile::<Words>::new(dir.join("none.txt"));
    assert!(matches!(
        absent.provide(),
        Err(TableError::Unreadable { .. })
    ));
}

#[test]
fn a_cache_without_the_file_is_absent_and_a_malformed_file_is_refused() {
    let dir = directory("cache");
    let _ = std::fs::write(dir.join("bad.txt"), "bad");
    let missing = Cache::<Words>::new(&dir, "nothing.txt");
    assert!(matches!(missing.provide(), Err(TableError::Absent { .. })));
    let bad = Cache::<Words>::new(&dir, "bad.txt");
    assert!(matches!(bad.provide(), Err(TableError::Malformed { .. })));
}

#[test]
fn an_in_memory_table_is_shared_not_copied() {
    let memory = InMemory::new(Provided {
        identity: TableIdentity::supplied("w", "1", "test", "CC0-1.0"),
        table: Words(vec!["a".into()]),
    });
    let (Ok(borrowed), Ok(shared)) = (memory.provide(), memory.shared()) else {
        panic!("memory provides");
    };
    assert!(core::ptr::eq(borrowed, Arc::as_ptr(&shared)));
    assert_eq!(memory.source(), Source::Memory);
}

#[test]
fn errors_map_to_their_codes_and_missing_names_every_way_to_supply() {
    let missing = TableError::Missing {
        kind: "words",
        wanted: "nothing selected".into(),
        option: "--words",
    };
    let text = missing.to_string();
    assert!(text.contains("--words"));
    assert!(text.contains("molframe data fetch <entry> --manifest <file>"));
    let code = |error: &TableError| molframe_core::Diagnostic::from(error).code();
    assert_eq!(code(&missing), molframe_core::Code::E7102);
    let unreadable = TableError::Unreadable {
        path: "p".into(),
        reason: "r".into(),
    };
    assert_eq!(code(&unreadable), molframe_core::Code::E7101);
    let malformed = TableError::Malformed {
        path: "p".into(),
        error: TableTextError::at("words table", 2, "bad"),
    };
    assert_eq!(code(&malformed), molframe_core::Code::E2201);
}
