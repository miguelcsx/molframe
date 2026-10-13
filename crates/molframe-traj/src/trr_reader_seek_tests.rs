//! Index scan, seek and sharing for the TRR reader.
#![allow(clippy::cast_precision_loss)]

use super::*;
use crate::{FrameIndex, TrrWriteOptions, write_trr};

/// Frames alternate between position-only and position-plus-velocity records.
fn fixture(frames: usize) -> (tempfile::TempDir, std::path::PathBuf) {
    let steps: Vec<Timestep> = (0..frames)
        .map(|i| Timestep {
            frame: i,
            time: Some(i as f64 * 0.25),
            positions: vec![[i as f32, 2.0, 3.0]; 9],
            velocities: (i % 3 == 0).then(|| vec![[0.5; 3]; 9]),
            ..Timestep::default()
        })
        .collect();
    let directory = tempfile::tempdir().expect("directory");
    let path = directory.path().join("s.trr");
    std::fs::write(
        &path,
        write_trr(&steps, TrrWriteOptions::default()).expect("encode"),
    )
    .expect("write");
    (directory, path)
}

fn sequential(path: &std::path::Path) -> Vec<Timestep> {
    let mut reader = TrrReader::open(path, 1 << 20).expect("open");
    let mut frame = Timestep::default();
    let mut out = Vec::new();
    while reader.read_next(&mut frame).expect("read") {
        out.push(frame.clone());
    }
    out
}

#[test]
fn seek_matches_sequential_for_variable_frames() {
    let (_directory, path) = fixture(40);
    let expected = sequential(&path);
    let mut reader = TrrReader::open(&path, 1 << 20).expect("open");
    assert_eq!(reader.random_access(), RandomAccess::ViaIndex);
    assert_eq!(reader.n_frames(), None);
    let mut frame = Timestep::default();
    for target in [39, 0, 20, 21, 1, 38, 3] {
        reader.seek(target).expect("seek");
        assert_eq!(reader.random_access(), RandomAccess::Full);
        assert!(reader.read_next(&mut frame).expect("read"));
        assert_eq!(frame, expected[target], "frame {target}");
    }
    reader.seek(40).expect("seek end");
    assert!(!reader.read_next(&mut frame).expect("eof"));
    assert_eq!(reader.n_frames(), Some(40));
}

#[test]
fn seek_then_continue_reads_following_frames() {
    let (_directory, path) = fixture(12);
    let expected = sequential(&path);
    let mut reader = TrrReader::open(&path, 1 << 20).expect("open");
    reader.seek(5).expect("seek");
    let mut frame = Timestep::default();
    for want in &expected[5..] {
        assert!(reader.read_next(&mut frame).expect("read"));
        assert_eq!(&frame, want);
    }
    assert!(!reader.read_next(&mut frame).expect("eof"));
}

#[test]
fn shared_index_is_checked_and_adopted() {
    let (_directory, path) = fixture(10);
    let index = TrrReader::open(&path, 1 << 20)
        .expect("open")
        .frame_index()
        .expect("index");
    assert_eq!(index.len(), 10);
    let adopted = TrrReader::open(&path, 1 << 20)
        .expect("open")
        .with_frame_index(&index)
        .expect("adopt");
    assert_eq!(adopted.n_frames(), Some(10));
    let (_other, other_path) = fixture(11);
    assert!(
        TrrReader::open(&other_path, 1 << 20)
            .expect("open")
            .with_frame_index(&index)
            .is_err()
    );
    let wrong = FrameIndex::new(FrameIndex::XTC, 9, index.stamp(), vec![0]).expect("index");
    assert!(
        TrrReader::open(&path, 1 << 20)
            .expect("open")
            .with_frame_index(&wrong)
            .is_err()
    );
}

#[test]
fn index_scan_respects_the_memory_limit() {
    let (_directory, path) = fixture(5000);
    let mut reader = TrrReader::open(&path, BUFFER + size_of::<TrrReader>() + 1024).expect("open");
    assert!(matches!(
        reader.seek(10),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
}

#[test]
fn truncated_file_is_refused_by_the_scan() {
    let (directory, path) = fixture(6);
    let bytes = std::fs::read(&path).expect("read");
    let cut = directory.path().join("cut.trr");
    std::fs::write(&cut, &bytes[..bytes.len() - 10]).expect("write");
    let mut reader = TrrReader::open(&cut, 1 << 20).expect("open");
    assert!(reader.seek(2).is_err());
}
