use std::path::Path;

use pdbiox_core::structure::UnitCell;
use tempfile::tempdir;

use crate::Timestep;

use super::{TngCompression, TngWriteOptions, parse_tng, write_tng};

fn frame(step: usize, time: f64, shift: f32) -> Timestep {
    Timestep {
        frame: step,
        time: Some(time),
        dt: Some(1.0),
        positions: vec![[1.0 + shift, 2.0, 3.0], [4.0, 5.0 + shift, 6.0]],
        velocities: Some(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
        forces: Some(vec![[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]),
        cell: Some(UnitCell {
            lengths: [20.0, 21.0, 22.0],
            angles: [90.0, 90.0, 90.0],
        }),
        ..Timestep::default()
    }
}

#[test]
fn lossless_round_trip_preserves_all_standard_frame_blocks() {
    let directory = tempdir().expect("temporary directory");
    let path = directory.path().join("complete.tng");
    let frames = [frame(0, 0.0, 0.0), frame(5, 2.0, 1.0)];

    write_tng(&path, &frames, TngWriteOptions::default()).expect("write TNG");
    let parsed = parse_tng(&path).expect("read TNG");

    assert_eq!(parsed.steps, vec![0, 5]);
    assert_eq!(parsed.frames.len(), 2);
    assert_eq!(parsed.frames[0].positions, frames[0].positions);
    let velocities = parsed.frames[1].velocities.as_ref().expect("velocities");
    assert!((velocities[1][0] - 0.4).abs() < 1.0e-6);
    assert_eq!(parsed.frames[0].forces, frames[0].forces);
    let cell = parsed.frames[1].cell.expect("cell");
    assert!((cell.lengths[1] - 21.0).abs() < 1.0e-5);
    assert!((parsed.frames[1].time.expect("time") - 2.0).abs() < 1.0e-9);
}

#[test]
fn angstrom_files_are_converted_at_the_boundary() {
    let directory = tempdir().expect("temporary directory");
    let path = directory.path().join("nanometres.tng");
    let frames = [frame(0, 0.0, 0.0)];
    let options = TngWriteOptions {
        distance_unit_exponent: -10,
        ..TngWriteOptions::default()
    };

    write_tng(&path, &frames, options).expect("write TNG");
    let parsed = parse_tng(&path).expect("read TNG");

    assert_eq!(parsed.distance_unit_exponent, -10);
    assert!((parsed.frames[0].positions[0][0] - 1.0).abs() < 1.0e-5);
    let forces = parsed.frames[0].forces.as_ref().expect("forces");
    assert!((forces[1][2] - 60.0).abs() < 1.0e-4);
}

#[test]
fn native_lossy_compression_is_readable() {
    let directory = tempdir().expect("temporary directory");
    let path = directory.path().join("lossy.tng");
    let frames = [frame(0, 0.0, 0.0), frame(1, 1.0, 0.25)];
    let options = TngWriteOptions {
        compression: TngCompression::Lossy { precision: 1_000.0 },
        ..TngWriteOptions::default()
    };

    write_tng(&path, &frames, options).expect("write compressed TNG");
    let parsed = parse_tng(&path).expect("read compressed TNG");

    assert_eq!(parsed.frames.len(), frames.len());
    assert!((parsed.frames[1].positions[0][0] - 1.25).abs() < 0.002);
}

#[test]
fn changing_atom_counts_are_rejected() {
    let directory = tempdir().expect("temporary directory");
    let path = directory.path().join("changing.tng");
    let mut frames = vec![frame(0, 0.0, 0.0), frame(1, 1.0, 0.0)];
    frames[1].positions.pop();

    assert!(write_tng(&path, &frames, TngWriteOptions::default()).is_err());
}

#[test]
fn external_corpus_reads_positions_cells_and_optional_vectors() {
    let Ok(root) = std::env::var("PDBIOX_TNG_CORPUS") else {
        return;
    };
    let plain = parse_tng(&Path::new(&root).join("argon_npt_compressed.tng"))
        .expect("read compressed corpus");
    let complete = parse_tng(&Path::new(&root).join("argon_npt_compressed_vels_forces.tng"))
        .expect("read complete corpus");

    assert_eq!(plain.frames.len(), 101);
    assert_eq!(plain.frames[0].positions.len(), 1_000);
    assert!(plain.frames[0].cell.is_some());
    assert!(plain.frames[0].velocities.is_none());
    assert_eq!(complete.frames.len(), 51);
    assert!(complete.frames[0].velocities.is_some());
    assert!(complete.frames[0].forces.is_some());
}
