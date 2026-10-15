use super::*;
use crate::{TrrWriteOptions, parse_trr, write_trr};

#[test]
fn pull_trr_matches_materialized_positions_and_auxiliary_arrays() {
    for precision in [TrrPrecision::Single, TrrPrecision::Double] {
        let frames: Vec<_> = (0_u16..3)
            .map(|i| Timestep {
                frame: usize::from(i),
                time: Some(f64::from(i) * 0.5),
                positions: vec![[1.0 + f32::from(i), 2.0, 3.0]; 12],
                velocities: Some(vec![[0.1; 3]; 12]),
                forces: Some(vec![[2.0; 3]; 12]),
                ..Timestep::default()
            })
            .collect();
        let bytes = write_trr(&frames, TrrWriteOptions { precision }).expect("encode");
        let expected = parse_trr(&bytes).expect("parse");
        let directory = tempfile::tempdir().expect("directory");
        let path = directory.path().join("sample.trr");
        std::fs::write(&path, bytes).expect("write");
        let mut reader = TrrReader::open(&path, 65536).expect("open");
        let mut frame = Timestep::default();
        let mut pointer = None;
        for expected in &expected.frames {
            assert!(reader.read_next(&mut frame).expect("frame"));
            assert_eq!(frame.positions, expected.positions);
            assert_eq!(frame.velocities, expected.velocities);
            assert_eq!(frame.forces, expected.forces);
            assert_eq!(frame.time, expected.time);
            assert_eq!(frame.dt, expected.dt);
            if let Some(pointer) = pointer {
                assert_eq!(frame.positions.as_ptr(), pointer);
            }
            pointer = Some(frame.positions.as_ptr());
        }
        assert!(!reader.read_next(&mut frame).expect("eof"));
    }
}

#[test]
fn insufficient_frame_budget_can_be_retried_without_losing_a_frame() {
    let frame = Timestep {
        positions: vec![[1.0; 3]; 100],
        ..Timestep::default()
    };
    let bytes = write_trr(&[frame], TrrWriteOptions::default()).expect("encode");
    let directory = tempfile::tempdir().expect("directory");
    let path = directory.path().join("sample.trr");
    std::fs::write(&path, bytes).expect("write");
    let mut reader = TrrReader::open(&path, 65536).expect("open");
    let mut frame = Timestep::default();
    assert!(matches!(
        reader.read_next_bounded(&mut frame, 9000),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
    assert!(frame.positions.is_empty());
    assert!(reader.read_next_bounded(&mut frame, 65536).expect("retry"));
    assert_eq!(frame.frame, 0);
    assert_eq!(frame.positions.len(), 100);
}

#[test]
fn pull_trr_decodes_frames_longer_than_one_read_block_in_both_precisions() {
    for precision in [TrrPrecision::Single, TrrPrecision::Double] {
        let atoms = 2 * BLOCK_ATOMS + 37;
        let ramp = |scale: f32| -> Vec<[f32; 3]> {
            (0_u16..)
                .take(atoms)
                .map(|i| {
                    let i = f32::from(i);
                    [i * scale, -i * scale, 0.5 * i * scale]
                })
                .collect()
        };
        let frames: Vec<_> = (0_u16..2)
            .map(|i| Timestep {
                frame: usize::from(i),
                time: Some(f64::from(i)),
                positions: ramp(0.25 + f32::from(i)),
                velocities: Some(ramp(0.125)),
                forces: Some(ramp(2.0)),
                ..Timestep::default()
            })
            .collect();
        let bytes = write_trr(&frames, TrrWriteOptions { precision }).expect("encode");
        let expected = parse_trr(&bytes).expect("parse");
        let directory = tempfile::tempdir().expect("directory");
        let path = directory.path().join("long.trr");
        std::fs::write(&path, &bytes).expect("write");
        let mut reader = TrrReader::open(&path, 1 << 20).expect("open");
        let mut frame = Timestep::default();
        for expected in &expected.frames {
            assert!(reader.read_next(&mut frame).expect("frame"));
            assert_eq!(frame.positions, expected.positions);
            assert_eq!(frame.velocities, expected.velocities);
            assert_eq!(frame.forces, expected.forces);
        }
        assert!(!reader.read_next(&mut frame).expect("eof"));
    }
}

#[test]
fn pull_trr_refuses_a_last_frame_cut_inside_its_positions() {
    let frames: Vec<_> = (0_u16..2)
        .map(|i| Timestep {
            frame: usize::from(i),
            time: Some(f64::from(i)),
            positions: vec![[1.0, 2.0, 3.0]; 300],
            ..Timestep::default()
        })
        .collect();
    let bytes = write_trr(&frames, TrrWriteOptions::default()).expect("encode");
    let directory = tempfile::tempdir().expect("directory");
    let path = directory.path().join("cut.trr");
    std::fs::write(&path, &bytes[..bytes.len() - 100]).expect("write");
    let mut reader = TrrReader::open(&path, 1 << 20).expect("open");
    let mut frame = Timestep::default();
    assert!(reader.read_next(&mut frame).expect("whole first frame"));
    assert!(reader.read_next(&mut frame).is_err());
}
