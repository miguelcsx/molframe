use std::path::{Path, PathBuf};

use super::XtcReader;
use crate::{
    RandomAccess, Timestep, TrajectoryError, TrajectoryFormat, TrajectoryIoError, TrajectoryReader,
    TrajectoryReaderOptions, XtcError, XtcWriteOptions, write_xtc,
};

fn fixture(frame_count: usize, atom_count: usize) -> (tempfile::TempDir, PathBuf) {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("stream.xtc");
    let frames = (0..frame_count)
        .map(|frame| {
            let frame_value = u16::try_from(frame).expect("fixture frame fits u16");
            Timestep {
                frame,
                time: Some(f64::from(frame_value) * 0.5),
                positions: (0..atom_count)
                    .map(|atom| {
                        let atom_value = u16::try_from(atom).expect("fixture atom fits u16");
                        [f32::from(frame_value), f32::from(atom_value), 1.0]
                    })
                    .collect(),
                ..Timestep::default()
            }
        })
        .collect::<Vec<_>>();
    let bytes = write_xtc(&frames, XtcWriteOptions::default()).expect("encode XTC fixture");
    std::fs::write(&path, bytes).expect("write XTC fixture");
    (directory, path)
}

fn open(path: &Path) -> XtcReader {
    XtcReader::open(path).expect("open XTC reader")
}

#[test]
fn ten_thousand_frames_reuse_one_output_and_one_workspace() {
    let (_directory, path) = fixture(10_000, 12);
    let mut reader = open(&path);
    assert_eq!(reader.n_frames(), None);
    let mut frame = Timestep::default();
    let mut pointer = None;
    let mut capacity = None;
    let mut workspace = None;
    let mut count = 0;
    while reader.read_next_frame(&mut frame).expect("read XTC frame") {
        let current_pointer = frame.positions.as_ptr();
        let current_capacity = frame.positions.capacity();
        if let Some(first) = pointer {
            assert_eq!(current_pointer, first);
        } else {
            pointer = Some(current_pointer);
        }
        if let Some(first) = capacity {
            assert_eq!(current_capacity, first);
        } else {
            capacity = Some(current_capacity);
        }
        let current_workspace = reader.workspace_bytes();
        if let Some(first) = workspace {
            assert_eq!(current_workspace, first);
        } else {
            workspace = Some(current_workspace);
        }
        count += 1;
    }
    assert_eq!(count, 10_000);
    assert_eq!(reader.n_frames(), Some(10_000));
    assert_eq!(reader.random_access(), RandomAccess::ViaIndex);
    assert!(reader.workspace_bytes() < 100_000);
}

#[test]
fn boxed_dispatch_is_object_safe_and_builds_offsets_only_on_seek() {
    let (_directory, path) = fixture(64, 12);
    let mut reader = crate::read_trajectory(&path, &TrajectoryReaderOptions::default())
        .expect("dispatch XTC pull reader");
    assert_eq!(reader.format(), "xtc");
    assert_eq!(reader.n_frames(), None);
    assert_eq!(reader.random_access(), RandomAccess::ViaIndex);
    reader.seek(42).expect("build bounded index and seek");
    assert_eq!(reader.n_frames(), Some(64));
    assert_eq!(reader.random_access(), RandomAccess::Full);
    let mut frame = Timestep::default();
    assert!(reader.read_next(&mut frame).expect("read indexed frame"));
    assert_eq!(frame.frame, 42);
    assert!((frame.positions[0][0] - 42.0).abs() < 0.01);
    assert_eq!(frame.dt, None);
}

#[test]
fn index_over_budget_degrades_to_forward_stream_without_losing_position() {
    let (_directory, path) = fixture(512, 1);
    let mut reader =
        XtcReader::open_with_memory_limit(&path, 68_000).expect("open tight XTC reader");
    let mut frame = Timestep::default();
    assert!(reader.read_next(&mut frame).expect("read first frame"));
    assert!(matches!(
        reader.seek(400),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
    assert_eq!(reader.random_access(), RandomAccess::None);
    assert_eq!(
        reader.seek(0),
        Err(TrajectoryError::RandomAccessUnavailable)
    );
    assert!(
        reader
            .read_next_frame(&mut frame)
            .expect("continue sequential read")
    );
    assert_eq!(frame.frame, 1);
}

#[test]
fn the_memory_contract_has_a_100_mb_default_and_no_maximum() {
    let (_directory, path) = fixture(1, 1);
    let _reader = open(&path);
    assert_eq!(
        TrajectoryReaderOptions::DEFAULT_MEMORY_LIMIT_BYTES,
        100_000_000
    );
    assert!(matches!(
        XtcReader::open_with_memory_limit(&path, 0),
        Err(XtcError::InvalidMemoryLimit { .. })
    ));
    assert!(
        XtcReader::open_with_memory_limit(&path, 64_000_000_000).is_ok(),
        "a caller who has provisioned the machine sets the ceiling, not the library"
    );
}

#[test]
fn pull_dispatch_refuses_to_promise_an_unimplemented_reader() {
    let result = crate::read_trajectory(
        Path::new("trajectory.tng"),
        &TrajectoryReaderOptions {
            format: Some(TrajectoryFormat::Tng),
            ..TrajectoryReaderOptions::default()
        },
    );
    assert!(matches!(
        result,
        Err(TrajectoryIoError::PullReaderUnavailable {
            format: TrajectoryFormat::Tng
        })
    ));
}
