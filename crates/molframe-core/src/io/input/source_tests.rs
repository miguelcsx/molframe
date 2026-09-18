use super::*;
use crate::{MemoryBudget, ScratchPolicy};
use std::io::{Seek, SeekFrom, Write};

#[test]
fn resident_windows_preserve_global_offsets() {
    let mut source = InputBuffer::from_bytes(b"abcdefgh".to_vec());
    let window = source.window(3, 3).expect("resident window");
    assert_eq!(window.start(), 3);
    assert_eq!(window.end(), Ok(6));
    assert_eq!(window.bytes(), b"def");
    assert_eq!(source.len_hint(), Some(8));
}

#[test]
fn window_end_rejects_u64_overflow_without_saturation() {
    let window = ByteWindow::new(u64::MAX, &[1]);
    assert_eq!(
        window.end().err().map(|error| error.code()),
        Some(Code::E1903)
    );
}

#[test]
fn a_file_window_is_bounded_by_the_shared_budget() {
    let budget = MemoryBudget::new(8).expect("non-zero budget");
    let context = ExecutionContext::builder()
        .memory_budget(budget)
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let path = std::env::temp_dir().join(format!("molframe-window-{}", std::process::id()));
    std::fs::write(&path, b"0123456789").expect("fixture");
    let mut source = WindowedFile::open(&path, 8, &context).expect("windowed file");
    assert_eq!(source.capacity(), 8);
    assert_eq!(source.window(4, 4).expect("window").bytes(), b"4567");
    assert!(source.window(0, 9).is_err());
    std::fs::remove_file(path).expect("fixture cleanup");
}

#[test]
fn sparse_offsets_past_u32_are_exact() {
    // The subject is 64-bit offset arithmetic, which any length past `u32::MAX`
    // exercises. The nominal 16 TiB is what a whole-file test would want, but
    // ext4 caps a file at exactly that, so the fixture takes the largest length
    // the filesystem actually backs.
    const SIXTEEN_TIB: u64 = 16 * 1024 * 1024 * 1024 * 1024;
    let context = ExecutionContext::default();
    let path = std::env::temp_dir().join(format!("molframe-sparse-window-{}", std::process::id()));
    let mut fixture = std::fs::File::create(&path).expect("sparse fixture");
    let Some(length) = largest_sparse_length(&fixture, SIXTEEN_TIB) else {
        println!("skipped: this filesystem cannot back a sparse file past u32::MAX");
        std::fs::remove_file(path).expect("sparse fixture cleanup");
        return;
    };
    fixture
        .seek(SeekFrom::Start(length - 1))
        .expect("sparse tail seek");
    fixture.write_all(b"x").expect("sparse tail write");
    drop(fixture);

    let mut source = WindowedFile::open(&path, 4096, &context).expect("windowed sparse file");
    assert_eq!(source.len_hint(), Some(length));
    let tail = source.window(length - 1, 1).expect("sparse tail window");
    assert_eq!(tail.start(), length - 1);
    assert_eq!(tail.end(), Ok(length));
    assert_eq!(tail.bytes(), b"x");
    std::fs::remove_file(path).expect("sparse fixture cleanup");
}

/// Returns the largest length at or below `target` that `fixture` can take.
///
/// Halves `target` until the filesystem accepts it, and returns `None` when not
/// even a length past `u32::MAX` can be created — the one case where the
/// property under test cannot be exercised at all.
fn largest_sparse_length(fixture: &std::fs::File, target: u64) -> Option<u64> {
    let mut length = target;
    while length > u64::from(u32::MAX) {
        if fixture.set_len(length).is_ok() {
            return Some(length);
        }
        length /= 2;
    }
    None
}
