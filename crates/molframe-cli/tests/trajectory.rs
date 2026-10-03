//! End-to-end streaming trajectory analysis on files written with known motion.

mod common;
use common::rows;
use molframe::trajectory::{DcdWriteOptions, Timestep, write_dcd};
use std::path::{Path, PathBuf};

/// Four frames of two atoms: the first stays put, the second oscillates by
/// +-1 Angstrom along x, so its RMSF is exactly 1 and the first atom's is 0.
fn oscillating(directory: &Path) -> PathBuf {
    let frames: Vec<Timestep> = [-1.0_f32, 1.0, -1.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(frame, offset)| Timestep {
            frame,
            positions: vec![[5.0, 5.0, 5.0], [10.0 + offset, 0.0, 0.0]],
            ..Timestep::default()
        })
        .collect();
    let bytes = write_dcd(&frames, &DcdWriteOptions::default()).expect("encode dcd");
    let path = directory.join("oscillating.dcd");
    std::fs::write(&path, bytes).expect("write dcd");
    path
}

#[test]
fn rmsf_is_zero_for_a_fixed_atom_and_the_amplitude_for_an_oscillating_one() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let file = oscillating(directory.path()).display().to_string();
    let table = rows(directory.path(), &["trajectory", "rmsf", &file]);
    assert_eq!(table[0], ["atom", "rmsf_angstrom"]);
    assert_eq!(table.len(), 3);
    let fixed: f64 = table[1][1].parse().expect("rmsf");
    let moving: f64 = table[2][1].parse().expect("rmsf");
    assert!(fixed.abs() < 1e-9, "{fixed}");
    assert!((moving - 1.0).abs() < 1e-6, "{moving}");
}

#[test]
fn a_trajectory_in_another_format_is_refused() {
    let directory = tempfile::tempdir().expect("scratch directory");
    let other = directory.path().join("frames.gro");
    std::fs::write(&other, "").expect("write file");
    let output = common::run(
        directory.path(),
        &["trajectory", "rmsf", other.to_str().expect("path")],
    );
    assert_eq!(output.status.code(), Some(7), "refused");
}
