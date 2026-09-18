use super::*;
use crate::ExecutionContext;
use std::fs::OpenOptions;
use std::io::{Seek, SeekFrom, Write};
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_root(label: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock must follow the Unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("molframe-{label}-{}-{nanos}", std::process::id()))
}

fn context(root: &std::path::Path, bytes: u64) -> ExecutionContext {
    ExecutionContext::builder()
        .temp_storage_policy(TempStoragePolicy::directory(root, bytes))
        .build()
        .expect("spill policy should be valid")
}

#[test]
fn spill_round_trip_reuses_one_payload_buffer() {
    let root = unique_root("spill-round-trip");
    let context = context(&root, 1024);
    let mut writer = context
        .create_spill_file("partition-0001")
        .expect("spill should be created");
    writer.append(b"alpha").expect("first record should fit");
    writer
        .append(b"longer-beta")
        .expect("second record should fit");
    let artifact = writer.finish().expect("spill should flush");
    assert_eq!(artifact.records(), 2);
    assert_eq!(context.spill_bytes(), artifact.bytes());

    let mut reader = artifact.reader().expect("artifact should open");
    let mut payload = Vec::with_capacity(32);
    assert_eq!(
        reader.read_next(&mut payload, 32).expect("record one"),
        Some(5)
    );
    assert_eq!(payload, b"alpha");
    assert_eq!(
        reader.read_next(&mut payload, 32).expect("record two"),
        Some(11)
    );
    assert_eq!(payload, b"longer-beta");
    assert_eq!(reader.read_next(&mut payload, 32).expect("clean eof"), None);

    let path = artifact.path().to_owned();
    drop(artifact);
    assert_eq!(context.spill_bytes(), 0);
    assert!(!path.exists());
    std::fs::remove_dir_all(root).expect("empty spill root should be removable");
}

#[test]
fn spill_budget_rejects_a_record_without_partial_bytes() {
    let root = unique_root("spill-budget");
    let context = context(&root, FILE_HEADER_BYTES + RECORD_HEADER_DISK_BYTES + 3);
    let mut writer = context
        .create_spill_file("bounded")
        .expect("header should fit");
    let error = writer
        .append(b"four")
        .expect_err("record should exceed disk ceiling");
    assert!(matches!(error, SpillError::BudgetExceeded { .. }));
    assert_eq!(context.spill_bytes(), FILE_HEADER_BYTES);
    drop(writer);
    assert_eq!(context.spill_bytes(), 0);
    std::fs::remove_dir_all(root).expect("empty spill root should be removable");
}

#[test]
fn spill_reader_detects_payload_corruption() {
    let root = unique_root("spill-corrupt");
    let context = context(&root, 1024);
    let mut writer = context
        .create_spill_file("corrupt")
        .expect("spill should be created");
    writer.append(b"payload").expect("record should fit");
    let artifact = writer.finish().expect("spill should flush");
    let mut file = OpenOptions::new()
        .write(true)
        .open(artifact.path())
        .expect("artifact should be writable for corruption test");
    file.seek(SeekFrom::End(-1))
        .expect("last byte should exist");
    file.write_all(&[0]).expect("last byte should be replaced");
    drop(file);

    let mut reader = artifact.reader().expect("header remains valid");
    let error = reader
        .read_next(&mut Vec::new(), 32)
        .expect_err("checksum mismatch must be visible");
    assert!(matches!(error, SpillError::CorruptRecord));
    drop(artifact);
    std::fs::remove_dir_all(root).expect("empty spill root should be removable");
}

#[test]
fn spill_reader_applies_the_memory_limit_before_growing_the_buffer() {
    let root = unique_root("spill-read-limit");
    let context = context(&root, 1024);
    let mut writer = context
        .create_spill_file("read-limit")
        .expect("spill should be created");
    writer.append(b"four").expect("record should fit on disk");
    let artifact = writer.finish().expect("spill should flush");
    let mut reader = artifact.reader().expect("artifact should open");
    let mut payload = Vec::new();
    let error = reader
        .read_next(&mut payload, 3)
        .expect_err("record must not exceed the caller memory ceiling");
    assert!(matches!(error, SpillError::RecordExceedsLimit { .. }));
    assert_eq!(payload.capacity(), 0);
    assert_eq!(
        reader
            .read_next(&mut payload, 4)
            .expect("the bounded record should remain retryable"),
        Some(4)
    );
    assert_eq!(payload, b"four");
    drop(artifact);
    std::fs::remove_dir_all(root).expect("empty spill root should be removable");
}

#[test]
fn spill_keys_cannot_escape_the_configured_root() {
    let root = unique_root("spill-key");
    let context = context(&root, 1024);
    let error = context
        .create_spill_file("../escape")
        .expect_err("path traversal must be rejected");
    assert!(matches!(error, SpillError::InvalidKey));
}
