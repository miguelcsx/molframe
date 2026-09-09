use super::*;
use crate::{ScratchPolicy, TempStoragePolicy};
use flate2::Compression as GzipLevel;
use flate2::write::GzEncoder;
use std::io::Write;

#[test]
fn compressed_source_requires_explicit_spill() {
    let path = unique_path("disabled.gz");
    write_gzip(&path, b"abcdefgh");
    let error = WindowedSourceFile::open(&path, 8, &ExecutionContext::default())
        .expect_err("default context must not spill");
    assert_eq!(error.code(), Code::E1902);
    std::fs::remove_file(path).expect("remove fixture");
}

#[test]
fn gzip_spool_is_bounded_seekable_and_deleted_exactly_once() {
    let base = unique_path("enabled");
    let path = base.with_extension("gz");
    let spill_root = base.with_extension("spill");
    write_gzip(&path, b"abcdefghijklmnopqrstuvwxyz");
    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .temp_storage_policy(TempStoragePolicy::directory(&spill_root, 4096))
        .build()
        .expect("context");
    let mut source = WindowedSourceFile::open(&path, 8, &context).expect("spooled source");
    assert!(source.spill_bytes() > 26);
    assert_eq!(source.len_hint(), Some(26));
    assert_eq!(
        source.window(6, 8).expect("cross-record").bytes(),
        b"ghijklmn"
    );
    assert_eq!(source.window(0, 8).expect("rollback").bytes(), b"abcdefgh");
    assert!(context.spill_bytes() > 0);
    drop(source);
    assert_eq!(context.spill_bytes(), 0);
    assert!(
        std::fs::read_dir(&spill_root)
            .expect("spill directory")
            .next()
            .is_none()
    );
    std::fs::remove_file(path).expect("remove fixture");
    std::fs::remove_dir(spill_root).expect("remove spill root");
}

#[test]
#[cfg(feature = "zstd")]
fn zstd_spool_uses_the_same_bounded_window_contract() {
    let base = unique_path("zstd-enabled");
    let path = base.with_extension("zst");
    let spill_root = base.with_extension("spill");
    let file = File::create(&path).expect("create fixture");
    let mut encoder = zstd::stream::write::Encoder::new(file, 1).expect("zstd encoder");
    encoder
        .write_all(b"abcdefghijklmnopqrstuvwxyz")
        .expect("encode fixture");
    encoder.finish().expect("finish fixture");
    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .temp_storage_policy(TempStoragePolicy::directory(&spill_root, 4096))
        .build()
        .expect("context");
    let mut source = WindowedSourceFile::open(&path, 8, &context).expect("spooled source");
    assert_eq!(source.len_hint(), Some(26));
    assert_eq!(
        source.window(14, 8).expect("cross-record").bytes(),
        b"opqrstuv"
    );
    drop(source);
    assert_eq!(context.spill_bytes(), 0);
    std::fs::remove_file(path).expect("remove fixture");
    std::fs::remove_dir(spill_root).expect("remove spill root");
}

fn write_gzip(path: &Path, bytes: &[u8]) {
    let file = File::create(path).expect("create fixture");
    let mut encoder = GzEncoder::new(file, GzipLevel::fast());
    encoder.write_all(bytes).expect("encode fixture");
    encoder.finish().expect("finish fixture");
}

fn unique_path(suffix: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "pdbiox-source-file-{}-{suffix}",
        std::process::id()
    ))
}
