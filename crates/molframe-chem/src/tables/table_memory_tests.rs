use super::*;
use crate::test_data::Words;

#[test]
fn profile_lines_skip_comments_blanks_and_empty_values() {
    let lines: Vec<_> =
        profile_lines("# note\n\n  words : lab \nempty:\nno colon\nreach: r\n").collect();
    assert_eq!(lines, vec![("words", "lab"), ("reach", "r")]);
}

#[test]
fn the_first_profile_line_for_a_kind_wins() {
    let mut tables = MemoryTables::new("test");
    tables.set_profile("words: first\nwords: second\n");
    assert_eq!(tables.profile_entry("words").as_deref(), Some("first"));
}

#[test]
fn setting_a_profile_replaces_the_one_before() {
    let mut tables = MemoryTables::new("test");
    tables.set_profile("words: first\n");
    tables.set_profile("other: x\n");
    assert_eq!(tables.profile_entry("words"), None);
}

#[test]
fn an_entry_is_found_by_name_and_extension_only() {
    let mut tables = MemoryTables::new("test");
    tables.insert("lab", "txt", TableBytes::from("a"));
    assert!(tables.entry("lab", "txt").is_some());
    assert!(tables.entry("lab", "reach").is_none());
    assert!(tables.entry("other", "txt").is_none());
}

#[test]
fn an_entry_decodes_once_per_kind() {
    let mut tables = MemoryTables::new("test");
    tables.insert("lab", "txt", TableBytes::from("a b"));
    let Some(entry) = tables.entry("lab", "txt") else {
        panic!("lab is held");
    };
    let (Ok(first), Ok(second)) = (entry.provided::<Words>(), entry.provided::<Words>()) else {
        panic!("lab decodes");
    };
    assert!(Arc::ptr_eq(&first, &second));
    assert_eq!(first.table, Words(vec!["a".into(), "b".into()]));
}

#[test]
fn a_malformed_entry_reports_the_same_error_each_time() {
    let mut tables = MemoryTables::new("test");
    tables.insert("lab", "txt", TableBytes::from("bad"));
    let Some(entry) = tables.entry("lab", "txt") else {
        panic!("lab is held");
    };
    let (Err(first), Err(second)) = (entry.provided::<Words>(), entry.provided::<Words>()) else {
        panic!("bad is refused");
    };
    assert!(matches!(first, TableError::Malformed { .. }));
    assert_eq!(first.to_string(), second.to_string());
    assert!(first.to_string().contains("test:lab.txt"));
}

#[test]
fn inserting_an_entry_again_replaces_it() {
    let mut tables = MemoryTables::new("test");
    tables.insert("lab", "txt", TableBytes::from("a"));
    tables.insert("lab", "txt", TableBytes::from("a b c"));
    let Some(Ok(provided)) = tables
        .entry("lab", "txt")
        .map(|entry| entry.provided::<Words>())
    else {
        panic!("lab decodes");
    };
    assert_eq!(provided.table.0.len(), 3);
}
