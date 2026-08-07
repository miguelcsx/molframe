use super::ExtensionStore;
use std::sync::Arc;

#[test]
fn values_require_both_matching_key_and_type() {
    let mut extensions = ExtensionStore::new();
    extensions.insert("example.count.v1", 7_u32);

    assert_eq!(extensions.get::<u32>("example.count.v1"), Some(&7));
    assert!(extensions.get::<u64>("example.count.v1").is_none());
    assert!(extensions.get::<u32>("example.other.v1").is_none());
}

#[test]
fn clones_share_values_and_keep_sorted_keys() {
    let mut extensions = ExtensionStore::new();
    extensions.insert("z.v1", String::from("shared"));
    extensions.insert("a.v1", 1_u8);
    let clone = extensions.clone();

    assert_eq!(clone.keys().collect::<Vec<_>>(), ["a.v1", "z.v1"]);
    let left = extensions.get_shared::<String>("z.v1");
    let right = clone.get_shared::<String>("z.v1");
    assert!(matches!((left, right), (Some(left), Some(right)) if Arc::ptr_eq(&left, &right)));
}

#[test]
fn clear_removes_every_value() {
    let mut extensions = ExtensionStore::new();
    extensions.insert("example.v1", true);
    assert_eq!(extensions.len(), 1);

    extensions.clear();

    assert!(extensions.is_empty());
}
