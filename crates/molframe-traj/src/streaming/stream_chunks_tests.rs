//! Parallel chunks must equal sequential chunks, in order, within bounds.

use super::*;
use crate::test_support::worker_context as context;
use crate::{DcdReader, DcdWriteOptions, XtcReader, XtcWriteOptions, write_dcd, write_xtc};
use molframe_core::CancellationToken;
use molframe_core::numeric::{usize_to_f32 as f32_of, usize_to_f64 as f64_of};
use std::path::{Path, PathBuf};

const FRAMES: usize = 700;

fn files() -> (tempfile::TempDir, PathBuf, PathBuf) {
    let steps: Vec<Timestep> = (0..FRAMES)
        .map(|i| Timestep {
            frame: i,
            time: Some(f64_of(i) * 0.5),
            positions: (0..6)
                .map(|a| [f32_of(i) * 0.25 + f32_of(a), f32_of(a), f32_of(i % 7)])
                .collect(),
            ..Timestep::default()
        })
        .collect();
    let dir = tempfile::tempdir().expect("dir");
    let dcd = dir.path().join("a.dcd");
    let xtc = dir.path().join("a.xtc");
    std::fs::write(
        &dcd,
        write_dcd(&steps, &DcdWriteOptions::default()).expect("dcd"),
    )
    .expect("write");
    std::fs::write(
        &xtc,
        write_xtc(&steps, XtcWriteOptions::default()).expect("xtc"),
    )
    .expect("write");
    (dir, dcd, xtc)
}

fn sequential<R: TrajectoryReader>(mut reader: R, range: Range<usize>) -> DenseChunk {
    reader.seek(range.start).expect("seek");
    read_dense_chunk(&mut reader, &mut Timestep::default(), range.len()).expect("chunk")
}

fn merged(chunks: impl Iterator<Item = Result<DenseChunk, TrajectoryError>>) -> DenseChunk {
    let mut out = DenseChunk::default();
    for chunk in chunks {
        let chunk = chunk.expect("chunk");
        out.atoms = chunk.atoms;
        out.indices.extend(chunk.indices);
        out.coordinates.extend(chunk.coordinates);
        out.times.extend(chunk.times);
    }
    out
}

fn same(left: &DenseChunk, right: &DenseChunk) {
    assert_eq!(left.indices, right.indices);
    assert_eq!(left.coordinates, right.coordinates);
    assert_eq!(left.times, right.times);
}

fn xtc(path: &Path) -> impl Fn() -> Result<XtcReader, TrajectoryError> + '_ {
    move || XtcReader::open(path).map_err(|_| TrajectoryError::RandomAccessUnavailable)
}

fn dcd(path: &Path) -> impl Fn() -> Result<DcdReader, TrajectoryError> + '_ {
    move || DcdReader::open(path).map_err(|_| TrajectoryError::RandomAccessUnavailable)
}

#[test]
fn every_worker_count_and_chunk_size_matches_the_sequential_read() {
    let (_dir, dcd_path, xtc_path) = files();
    let wanted_xtc = sequential(xtc(&xtc_path)().expect("open"), 0..FRAMES);
    let wanted_dcd = sequential(dcd(&dcd_path)().expect("open"), 0..FRAMES);
    for workers in [1, 2, 3, 8] {
        for chunk in [1, 7, 64, 1000] {
            let ctx = context(workers);
            let got = ParallelChunks::open(xtc(&xtc_path), 0..FRAMES, chunk, &ctx).expect("open");
            same(&merged(got), &wanted_xtc);
            let got = ParallelChunks::open(dcd(&dcd_path), 0..FRAMES, chunk, &ctx).expect("open");
            same(&merged(got), &wanted_dcd);
        }
    }
}

#[test]
fn chunks_are_whole_blocks_in_order_and_the_range_is_honoured() {
    let (_dir, _dcd, xtc_path) = files();
    let ctx = context(4);
    let chunks: Vec<DenseChunk> = ParallelChunks::open(xtc(&xtc_path), 100..350, 100, &ctx)
        .expect("open")
        .map(|chunk| chunk.expect("chunk"))
        .collect();
    let sizes: Vec<usize> = chunks.iter().map(DenseChunk::frames).collect();
    assert_eq!(sizes, [100, 100, 50]);
    assert_eq!(chunks[0].indices[0], 100);
    assert_eq!(chunks[2].indices[49], 349);
}

#[test]
fn ranges_past_the_end_are_clipped_or_empty() {
    let (_dir, _dcd, xtc_path) = files();
    let ctx = context(2);
    let tail = ParallelChunks::open(xtc(&xtc_path), 650..9000, 30, &ctx).expect("open");
    assert_eq!(merged(tail).frames(), 50);
    let none = ParallelChunks::open(xtc(&xtc_path), 900..950, 30, &ctx).expect("open");
    assert_eq!(none.count(), 0);
}

#[test]
fn dropping_early_stops_the_workers() {
    let (_dir, _dcd, xtc_path) = files();
    let ctx = context(4);
    let mut chunks = ParallelChunks::open(xtc(&xtc_path), 0..FRAMES, 5, &ctx).expect("open");
    assert!(chunks.next().expect("first").is_ok());
    drop(chunks);
    assert_eq!(ctx.reserved_bytes(), 0);
}

#[test]
fn cancellation_surfaces_as_an_error() {
    let (_dir, _dcd, xtc_path) = files();
    let token = CancellationToken::new();
    let ctx = context(2).with_cancellation(token.clone());
    token.cancel();
    let mut chunks = ParallelChunks::open(xtc(&xtc_path), 0..FRAMES, 10, &ctx).expect("open");
    assert!(matches!(
        chunks.next(),
        Some(Err(TrajectoryError::Cancelled))
    ));
    assert!(chunks.next().is_none());
}

#[test]
fn a_window_larger_than_the_budget_is_refused_before_any_decode() {
    let (_dir, _dcd, xtc_path) = files();
    let ctx = context(4);
    let _held = ctx.try_reserve(ctx.available_bytes() - 1000).expect("hold");
    assert!(matches!(
        ParallelChunks::open(xtc(&xtc_path), 0..FRAMES, 500, &ctx),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
}

#[test]
fn a_source_shorter_than_declared_is_an_error_not_a_hang() {
    let (dir, dcd_path, _xtc) = files();
    let bytes = std::fs::read(&dcd_path).expect("read");
    let cut = dir.path().join("cut.dcd");
    std::fs::write(&cut, &bytes[..bytes.len() - 40]).expect("write");
    let ctx = context(3);
    let outcome: Result<Vec<_>, _> = match ParallelChunks::open(dcd(&cut), 0..FRAMES, 50, &ctx) {
        Ok(chunks) => chunks.collect(),
        Err(error) => Err(error),
    };
    assert!(outcome.is_err());
}
