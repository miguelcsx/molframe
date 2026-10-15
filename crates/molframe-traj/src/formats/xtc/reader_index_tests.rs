//! Index scan, sharing and cache for the XTC reader.

use molframe_core::numeric::{usize_to_f32 as f32_of, usize_to_f64 as f64_of};
use std::path::PathBuf;

use super::XtcReader;
use crate::{
    FrameIndex, RandomAccess, Timestep, TrajectoryReader, XtcError, XtcWriteOptions, write_xtc,
};

/// Frames whose compressed size varies, so no stride can locate them.
fn fixture(name: &str, frames: usize) -> (tempfile::TempDir, PathBuf) {
    let directory = tempfile::tempdir().expect("directory");
    let path = directory.path().join(name);
    let steps: Vec<Timestep> = (0..frames)
        .map(|frame| Timestep {
            frame,
            time: Some(f64_of(frame) * 0.5),
            positions: (0..40)
                .map(|atom| {
                    let wild = f32_of((frame * 37 + atom * 11) % 17) * f32_of(frame % 5);
                    [f32_of(atom) * 0.1 + wild, f32_of(frame % 9), f32_of(atom)]
                })
                .collect(),
            ..Timestep::default()
        })
        .collect();
    std::fs::write(
        &path,
        write_xtc(&steps, XtcWriteOptions::default()).expect("encode"),
    )
    .expect("write");
    (directory, path)
}

fn sequential(path: &std::path::Path) -> Vec<Timestep> {
    let mut reader = XtcReader::open(path).expect("open");
    let mut frame = Timestep::default();
    let mut out = Vec::new();
    while reader.read_next(&mut frame).expect("read") {
        out.push(frame.clone());
    }
    out
}

#[test]
fn frame_sizes_vary_and_seek_matches_sequential_including_dt() {
    let (_directory, path) = fixture("a.xtc", 60);
    let expected = sequential(&path);
    let mut reader = XtcReader::open(&path).expect("open");
    let index = reader.frame_index().expect("index");
    let sizes: std::collections::BTreeSet<u64> =
        index.offsets().windows(2).map(|w| w[1] - w[0]).collect();
    assert!(sizes.len() > 1, "fixture must have variable frame sizes");
    let mut frame = Timestep::default();
    for target in [59, 0, 31, 30, 1, 58, 17] {
        reader.seek(target).expect("seek");
        assert!(reader.read_next(&mut frame).expect("read"));
        assert_eq!(frame, expected[target], "frame {target}");
    }
}

#[test]
fn shared_index_skips_the_scan_and_is_checked() {
    let (_directory, path) = fixture("a.xtc", 20);
    let index = XtcReader::open(&path)
        .expect("open")
        .frame_index()
        .expect("index");
    let adopted = XtcReader::open(&path)
        .expect("open")
        .with_frame_index(&index)
        .expect("adopt");
    assert_eq!(adopted.n_frames(), Some(20));
    assert_eq!(adopted.random_access(), RandomAccess::Full);
    let (_other_dir, other) = fixture("b.xtc", 21);
    assert!(matches!(
        XtcReader::open(&other)
            .expect("open")
            .with_frame_index(&index),
        Err(XtcError::Index(_))
    ));
    let forged = FrameIndex::new(FrameIndex::TRR, 40, index.stamp(), vec![0]).expect("forged");
    assert!(
        XtcReader::open(&path)
            .expect("open")
            .with_frame_index(&forged)
            .is_err()
    );
}

#[test]
fn cached_index_is_written_reused_and_rebuilt_when_stale() {
    let (directory, path) = fixture("a.xtc", 25);
    let cache = directory.path().join("a.xtc.idx");
    let first = XtcReader::open_cached(&path, 1 << 26, &cache).expect("open");
    assert_eq!(first.n_frames(), Some(25));
    let stored = FrameIndex::read(&cache).expect("cache written");
    assert_eq!(stored.len(), 25);

    let reused = XtcReader::open_cached(&path, 1 << 26, &cache).expect("open");
    assert_eq!(reused.n_frames(), Some(25));

    std::fs::write(&cache, b"not an index").expect("damage");
    let healed = XtcReader::open_cached(&path, 1 << 26, &cache).expect("open");
    assert_eq!(healed.n_frames(), Some(25));
    assert_eq!(FrameIndex::read(&cache).expect("rewritten"), stored);

    let (_dir, longer) = fixture("c.xtc", 30);
    let bytes = std::fs::read(&longer).expect("read");
    std::fs::write(&path, bytes).expect("replace source");
    let stale = XtcReader::open_cached(&path, 1 << 26, &cache).expect("open");
    assert_eq!(stale.n_frames(), Some(30));
}

#[test]
fn trait_sharing_adopts_only_a_matching_index() {
    let (_directory, path) = fixture("a.xtc", 20);
    let mut first = XtcReader::open(&path).expect("open");
    let index = first.share_index().expect("index");
    let mut second = XtcReader::open(&path).expect("open");
    assert!(second.adopt_index(&index));
    assert_eq!(second.n_frames(), Some(20));
    let (_other_dir, other) = fixture("b.xtc", 21);
    assert!(!XtcReader::open(&other).expect("open").adopt_index(&index));
}
