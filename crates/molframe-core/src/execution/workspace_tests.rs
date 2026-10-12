use super::*;
use crate::execution::TempStoragePolicy;

#[test]
fn native_tool_directories_share_the_spill_budget_and_release_it_after_cleanup() {
    let root = tempfile::tempdir().expect("root");
    let context = ExecutionContext::builder()
        .temp_storage_policy(TempStoragePolicy::directory(root.path(), 100))
        .build()
        .expect("context");
    let first = context.temporary_workspace(60).expect("first tool");
    let path = first.path().to_owned();
    assert_eq!(context.spill_bytes(), 60);
    assert!(context.clone().temporary_workspace(41).is_err());
    let second = context.temporary_workspace(40).expect("second tool");
    assert_ne!(first.path(), second.path());
    std::fs::write(first.path().join("output"), [0u8; 59]).expect("output");
    assert_eq!(first.check_size().expect("size"), 59);
    drop(first);
    assert!(!path.exists());
    assert_eq!(context.spill_bytes(), 40);
    drop(second);
    assert_eq!(context.spill_bytes(), 0);
}

#[test]
fn native_workspaces_reject_excess_output_and_do_not_enable_disabled_storage() {
    assert!(matches!(
        ExecutionContext::default().temporary_workspace(10),
        Err(SpillError::Disabled)
    ));
    let root = tempfile::tempdir().expect("root");
    let context = ExecutionContext::builder()
        .temp_storage_policy(TempStoragePolicy::directory(root.path(), 10))
        .build()
        .expect("context");
    let workspace = context.temporary_workspace(10).expect("workspace");
    std::fs::write(workspace.path().join("output"), [0u8; 11]).expect("output");
    assert!(matches!(
        workspace.check_size(),
        Err(SpillError::BudgetExceeded { .. })
    ));
    std::fs::remove_file(workspace.path().join("output")).expect("remove oversized output");
    std::fs::create_dir(workspace.path().join("nested")).expect("nested directory");
    assert!(matches!(workspace.check_size(), Err(SpillError::Io(_))));
}
