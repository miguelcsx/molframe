//! Seeking a DCD frame must decode exactly what a sequential read yields.
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap
)]

use super::DcdReader;
use crate::{DcdWriteOptions, Timestep, TrajectoryReader, write_dcd};
use molframe_core::structure::UnitCell;
use std::path::PathBuf;

fn written(frames: usize, cell: bool) -> (tempfile::TempDir, PathBuf, Vec<Timestep>) {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("seek.dcd");
    let steps: Vec<Timestep> = (0..frames)
        .map(|frame| Timestep {
            frame,
            time: Some(frame as f64),
            positions: (0..7)
                .map(|atom| [frame as f32, atom as f32, (frame * atom) as f32])
                .collect(),
            cell: cell.then_some(UnitCell {
                lengths: [10.0 + frame as f64, 11.0, 12.0],
                angles: [90.0; 3],
            }),
            ..Timestep::default()
        })
        .collect();
    let bytes = write_dcd(&steps, &DcdWriteOptions::default()).expect("encode");
    std::fs::write(&path, bytes).expect("write");
    (directory, path, steps)
}

fn sequential(path: &std::path::Path) -> Vec<Timestep> {
    let mut reader = DcdReader::open(path).expect("open");
    let mut out = Vec::new();
    let mut frame = Timestep::default();
    while reader.read_next(&mut frame).expect("read") {
        out.push(frame.clone());
    }
    out
}

#[test]
fn seeking_any_frame_matches_the_sequential_decode_with_and_without_cells() {
    for cell in [false, true] {
        let (_dir, path, _) = written(23, cell);
        let all = sequential(&path);
        let mut reader = DcdReader::open(&path).expect("open");
        let mut frame = Timestep::default();
        for target in [22, 0, 11, 5, 5, 17, 1] {
            reader.seek(target).expect("seek");
            assert!(reader.read_next(&mut frame).expect("read"));
            assert_eq!(frame, all[target], "frame {target} cell {cell}");
        }
        reader.seek(23).expect("seek to end");
        assert!(!reader.read_next(&mut frame).expect("eof"));
        assert!(reader.seek(24).is_err());
        reader.rewind().expect("rewind");
        assert!(reader.read_next(&mut frame).expect("read"));
        assert_eq!(frame.frame, 0);
    }
}

#[test]
fn seeking_a_fixed_atom_file_restores_the_frame_zero_atoms() {
    let (_dir, path) = fixed_file();
    let all = sequential(&path);
    assert_eq!(all.len(), 3);
    let mut reader = DcdReader::open(&path).expect("open");
    let mut frame = Timestep::default();
    for target in [2, 1, 0, 2] {
        reader.seek(target).expect("seek");
        assert!(reader.read_next(&mut frame).expect("read"));
        assert_eq!(frame, all[target]);
    }
}

fn record(out: &mut Vec<u8>, payload: &[u8]) {
    let length = (payload.len() as i32).to_le_bytes();
    out.extend_from_slice(&length);
    out.extend_from_slice(payload);
    out.extend_from_slice(&length);
}

fn floats(values: &[f32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

fn fixed_file() -> (tempfile::TempDir, PathBuf) {
    let mut bytes = Vec::new();
    let mut header = vec![0u8; 84];
    header[0..4].copy_from_slice(b"CORD");
    header[4..8].copy_from_slice(&3i32.to_le_bytes());
    header[12..16].copy_from_slice(&1i32.to_le_bytes());
    header[36..40].copy_from_slice(&1i32.to_le_bytes());
    header[40..44].copy_from_slice(&1.0f32.to_le_bytes());
    header[80..84].copy_from_slice(&24i32.to_le_bytes());
    record(&mut bytes, &header);
    let mut title = vec![b' '; 84];
    title[0..4].copy_from_slice(&1i32.to_le_bytes());
    record(&mut bytes, &title);
    record(&mut bytes, &3i32.to_le_bytes());
    record(
        &mut bytes,
        &[2i32.to_le_bytes(), 3i32.to_le_bytes()].concat(),
    );
    for axis in [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]] {
        record(&mut bytes, &floats(&axis));
    }
    for step in 1..3 {
        let shift = step as f32 * 10.0;
        for axis in [[20.0, 30.0], [50.0, 60.0], [80.0, 90.0]] {
            record(&mut bytes, &floats(&[axis[0] + shift, axis[1] + shift]));
        }
    }
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("fixed.dcd");
    std::fs::write(&path, bytes).expect("write");
    (directory, path)
}
