use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use super::DcdReader;
use crate::{
    DcdError, DcdWriteOptions, RandomAccess, Timestep, TrajectoryReader, TrajectoryReaderOptions,
    write_dcd,
};

fn fixture(frame_count: usize, atom_count: usize) -> (tempfile::TempDir, PathBuf) {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("stream.dcd");
    let frames = (0..frame_count)
        .map(|frame| Timestep {
            frame,
            time: Some(f64::from(
                u32::try_from(frame).expect("fixture frame fits in u32"),
            )),
            positions: (0..atom_count)
                .map(|atom| {
                    [
                        f32::from(u16::try_from(frame).expect("fixture frame fits in u16")),
                        f32::from(u16::try_from(atom).expect("fixture atom fits in u16")),
                        1.0,
                    ]
                })
                .collect(),
            ..Timestep::default()
        })
        .collect::<Vec<_>>();
    let bytes = write_dcd(&frames, &DcdWriteOptions::default()).expect("encode DCD fixture");
    std::fs::write(&path, bytes).expect("write DCD fixture");
    (directory, path)
}

fn open(path: &Path) -> DcdReader {
    DcdReader::open(path).expect("open DCD reader")
}

#[test]
fn ten_thousand_frames_reuse_two_coordinate_buffers() {
    let (_directory, path) = fixture(10_000, 12);
    let mut reader = open(&path);
    let mut frame = Timestep::default();
    let mut pointers = BTreeSet::new();
    let mut capacities = BTreeSet::new();
    let mut stable_workspace = None;
    let mut count = 0_u16;
    while reader.read_next_frame(&mut frame).expect("read DCD frame") {
        pointers.insert(frame.positions.as_ptr());
        capacities.insert(frame.positions.capacity());
        if count >= 1 {
            let current = reader.workspace_bytes();
            if let Some(expected) = stable_workspace {
                assert_eq!(current, expected);
            } else {
                stable_workspace = Some(current);
            }
        }
        assert!((frame.positions[0][0] - f32::from(count)).abs() <= f32::EPSILON);
        count += 1;
    }
    assert_eq!(count, 10_000);
    assert!(pointers.len() <= 2);
    assert_eq!(capacities, BTreeSet::from([12]));
    assert!(reader.workspace_bytes() < 100_000);
}

#[test]
fn workspace_is_independent_of_total_frame_count() {
    let (_short_directory, short_path) = fixture(1, 64);
    let (_long_directory, long_path) = fixture(2_048, 64);
    let short = open(&short_path);
    let long = open(&long_path);
    assert_eq!(short.workspace_bytes(), long.workspace_bytes());
    assert_eq!(short.n_frames(), Some(1));
    assert_eq!(long.n_frames(), Some(2_048));
}

#[test]
fn pull_dispatch_reads_dcd_without_materialising_all_frames() {
    let (_directory, path) = fixture(32, 8);
    let mut reader = crate::read_trajectory(&path, &TrajectoryReaderOptions::default())
        .expect("dispatch DCD pull reader");
    assert_eq!(reader.format(), "dcd");
    assert_eq!(reader.n_frames(), Some(32));
    assert_eq!(reader.random_access(), RandomAccess::None);
    let mut frame = Timestep::default();
    assert!(reader.read_next(&mut frame).expect("read dispatched DCD"));
    assert_eq!(frame.frame, 0);
    assert_eq!(frame.positions.len(), 8);
}

#[test]
fn memory_ceiling_rejects_a_frame_before_decoding() {
    let (_directory, path) = fixture(1, 1_024);
    let result = DcdReader::open_with_memory_limit(&path, 70_000);
    assert!(matches!(result, Err(DcdError::MemoryLimit { .. })));
}

#[test]
fn bounded_reader_matches_materialised_parser() {
    let (_directory, path) = fixture(17, 23);
    let bytes = std::fs::read(&path).expect("read fixture bytes");
    let expected = crate::parse_dcd(&bytes).expect("materialise DCD fixture");
    let mut reader = open(&path);
    let mut actual = Timestep::default();
    for expected_frame in expected.frames {
        assert!(reader.read_next_frame(&mut actual).expect("read DCD frame"));
        assert_eq!(actual, expected_frame);
    }
    assert!(!reader.read_next_frame(&mut actual).expect("reach DCD EOF"));
}

#[test]
fn fixed_atoms_are_retained_without_a_previous_frame_allocation() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("fixed.dcd");
    let endian = crate::DcdEndian::Little;
    let mut bytes = Vec::new();
    let mut header = vec![0_u8; 84];
    header[0..4].copy_from_slice(b"CORD");
    put_i32(&mut header, 4, 2);
    put_i32(&mut header, 12, 1);
    put_i32(&mut header, 36, 1);
    header[40..44].copy_from_slice(&1.0_f32.to_le_bytes());
    put_i32(&mut header, 80, 24);
    add_record(&mut bytes, &header, endian);
    let mut title = vec![b' '; 84];
    put_i32(&mut title, 0, 1);
    add_record(&mut bytes, &title, endian);
    add_record(&mut bytes, &3_i32.to_le_bytes(), endian);
    add_record(
        &mut bytes,
        &[2_i32.to_le_bytes(), 3_i32.to_le_bytes()].concat(),
        endian,
    );
    for axis in [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]] {
        add_record(&mut bytes, &floats(&axis), endian);
    }
    for axis in [[20.0, 30.0], [50.0, 60.0], [80.0, 90.0]] {
        add_record(&mut bytes, &floats(&axis), endian);
    }
    std::fs::write(&path, bytes).expect("write fixed DCD fixture");
    let mut reader = open(&path);
    let mut frame = Timestep::default();
    assert!(
        reader
            .read_next_frame(&mut frame)
            .expect("read first frame")
    );
    assert!(
        reader
            .read_next_frame(&mut frame)
            .expect("read second frame")
    );
    assert_position_close(frame.positions[0], [1.0, 4.0, 7.0]);
    assert_position_close(frame.positions[2], [30.0, 60.0, 90.0]);
}

fn assert_position_close(actual: [f32; 3], expected: [f32; 3]) {
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (*actual - expected).abs() <= f32::EPSILON)
    );
}

fn put_i32(bytes: &mut [u8], offset: usize, value: i32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn add_record(output: &mut Vec<u8>, payload: &[u8], endian: crate::DcdEndian) {
    let length = i32::try_from(payload.len()).expect("test record fits DCD marker");
    let marker = match endian {
        crate::DcdEndian::Little => length.to_le_bytes(),
        crate::DcdEndian::Big => length.to_be_bytes(),
    };
    output.extend(marker);
    output.extend(payload);
    output.extend(marker);
}

fn floats(values: &[f32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}
