//! Parallel range decoding must match sequential reads and respect its bounds.

use super::*;
use crate::shared_tests::worker_context as context;
use crate::{DcdReader, DcdWriteOptions, XtcReader, XtcWriteOptions, write_dcd, write_xtc};
use molframe_core::numeric::{usize_to_f32 as f32_of, usize_to_f64 as f64_of};
use std::path::{Path, PathBuf};

const FRAMES: usize = 1500;

fn frames() -> Vec<Timestep> {
    (0..FRAMES)
        .map(|i| Timestep {
            frame: i,
            time: Some(f64_of(i)),
            positions: (0..5)
                .map(|a| [f32_of(i) * 0.25, f32_of(a), f32_of(i % 7)])
                .collect(),
            ..Timestep::default()
        })
        .collect()
}

fn files() -> (tempfile::TempDir, PathBuf, PathBuf) {
    let dir = tempfile::tempdir().expect("dir");
    let dcd = dir.path().join("a.dcd");
    let xtc = dir.path().join("a.xtc");
    let steps = frames();
    std::fs::write(
        &dcd,
        write_dcd(&steps, &DcdWriteOptions::default()).expect("dcd"),
    )
    .expect("w");
    std::fs::write(
        &xtc,
        write_xtc(&steps, XtcWriteOptions::default()).expect("xtc"),
    )
    .expect("w");
    (dir, dcd, xtc)
}

fn collect<R: TrajectoryReader>(
    open: impl Fn() -> Result<R, TrajectoryError> + Sync,
    range: Range<usize>,
    workers: usize,
) -> Vec<(usize, f32)> {
    let context = context(workers);
    let mut out = Vec::new();
    let count = run_frame_ranges::<_, _, TrajectoryError, _, _, _>(
        open,
        range,
        &context,
        1 << 20,
        |frame, _| Ok(frame.positions[0][0] + frame.positions[4][2]),
        |ordinal, value| {
            out.push((ordinal, value));
            Ok(())
        },
    )
    .expect("ranges");
    assert_eq!(usize::try_from(count).expect("small"), out.len());
    assert_eq!(context.reserved_bytes(), 0);
    out
}

fn dcd(path: &Path) -> impl Fn() -> Result<DcdReader, TrajectoryError> + Sync + '_ {
    move || DcdReader::open(path).map_err(|_| TrajectoryError::RandomAccessUnavailable)
}

fn xtc(path: &Path) -> impl Fn() -> Result<XtcReader, TrajectoryError> + Sync + '_ {
    move || XtcReader::open(path).map_err(|_| TrajectoryError::RandomAccessUnavailable)
}

#[test]
fn every_worker_count_yields_the_serial_values_in_order() {
    let (_dir, dcd_path, xtc_path) = files();
    for range in [0..FRAMES, 300..1301, 1499..1500, 10..10] {
        let serial = collect(dcd(&dcd_path), range.clone(), 1);
        assert_eq!(serial.len(), range.len());
        for (position, (ordinal, _)) in serial.iter().enumerate() {
            assert_eq!(*ordinal, range.start + position);
        }
        for workers in [2, 3, 8] {
            assert_eq!(collect(dcd(&dcd_path), range.clone(), workers), serial);
            assert_eq!(
                collect(xtc(&xtc_path), range.clone(), workers).len(),
                serial.len()
            );
        }
    }
}

#[test]
fn xtc_parallel_decode_matches_serial_decode() {
    let (_dir, _, xtc_path) = files();
    let serial = collect(xtc(&xtc_path), 0..FRAMES, 1);
    assert_eq!(collect(xtc(&xtc_path), 0..FRAMES, 4), serial);
}

#[test]
fn the_range_is_clipped_to_the_frames_present() {
    let (_dir, dcd_path, _) = files();
    let got = collect(dcd(&dcd_path), 1400..usize::MAX, 4);
    assert_eq!(got.len(), 100);
    assert!(collect(dcd(&dcd_path), 5000..6000, 4).is_empty());
}

#[test]
fn a_compute_failure_emits_earlier_frames_then_reports_the_error() {
    let (_dir, dcd_path, _) = files();
    let context = context(4);
    let mut seen = Vec::new();
    let result = run_frame_ranges::<_, _, TrajectoryError, _, _, _>(
        dcd(&dcd_path),
        0..FRAMES,
        &context,
        1 << 20,
        |frame, _| {
            if frame.positions[0][0] > 100.0 {
                Err(TrajectoryError::Cancelled)
            } else {
                Ok(())
            }
        },
        |ordinal, ()| {
            seen.push(ordinal);
            Ok(())
        },
    );
    assert!(matches!(result, Err(TrajectoryError::Cancelled)));
    assert_eq!(seen, (0..seen.len()).collect::<Vec<_>>());
    assert!(seen.len() <= 401);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn a_memory_budget_too_small_for_the_workers_is_refused() {
    let (_dir, dcd_path, _) = files();
    let context = ExecutionContext::builder()
        .worker_budget(4)
        .memory_budget(molframe_core::MemoryBudget::new(16 << 20).expect("budget"))
        .build()
        .expect("context");
    let result = run_frame_ranges::<_, _, TrajectoryError, _, _, _>(
        dcd(&dcd_path),
        0..FRAMES,
        &context,
        6 << 20,
        |_, _| Ok(()),
        |_, ()| Ok(()),
    );
    assert!(matches!(result, Err(TrajectoryError::MemoryLimit { .. })));
}
