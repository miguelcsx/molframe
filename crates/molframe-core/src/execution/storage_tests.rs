use super::*;

#[test]
fn temporary_storage_is_disabled_by_default() {
    let policy = TempStoragePolicy::default();
    assert_eq!(policy.root(), None);
    assert_eq!(policy.max_bytes(), 0);
}

#[test]
fn directory_policy_retains_its_explicit_limits() {
    let policy = TempStoragePolicy::directory("spill", 42);
    assert_eq!(policy.root(), Some(Path::new("spill")));
    assert_eq!(policy.max_bytes(), 42);
}
