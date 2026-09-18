use super::*;
use molframe_core::{MemoryBudget, ScratchPolicy};

fn context(bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(bytes).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context")
}

#[test]
fn admission_failure_precedes_opening_the_input() {
    let context = context(100);
    let result = read_trajectory_in(
        Path::new("absent.trr"),
        &TrajectoryReaderOptions {
            memory_limit_bytes: 101,
            ..TrajectoryReaderOptions::default()
        },
        &context,
    );
    assert!(matches!(
        result,
        Err(TrajectoryIoError::Reader(
            TrajectoryError::MemoryLimit { .. }
        ))
    ));
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn an_open_error_releases_the_reader_allowance() {
    let context = context(32768);
    let directory = tempfile::tempdir().expect("directory");
    let result = read_trajectory_in(
        &directory.path().join("absent.trr"),
        &TrajectoryReaderOptions {
            memory_limit_bytes: 32768,
            ..TrajectoryReaderOptions::default()
        },
        &context,
    );
    assert!(result.is_err());
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn reader_and_frames_share_the_context_until_their_owners_release_them() {
    let context = context(65536);
    let directory = tempfile::tempdir().expect("directory");
    let path = directory.path().join("frames.trr");
    let frames = [Timestep {
        positions: vec![[1.0; 3]; 10],
        ..Timestep::default()
    }];
    std::fs::write(
        &path,
        crate::write_trr(&frames, crate::TrrWriteOptions::default()).expect("encode"),
    )
    .expect("write");
    let mut reader = read_trajectory_in(
        &path,
        &TrajectoryReaderOptions {
            memory_limit_bytes: 32768,
            ..TrajectoryReaderOptions::default()
        },
        &context,
    )
    .expect("open");
    assert_eq!(context.reserved_bytes(), 32768);
    let count = crate::run_analysis_stream(&mut reader, &context, 32768, |frame, context| {
        assert_eq!(context.reserved_bytes(), 65536);
        assert_eq!(frame.positions, frames[0].positions);
        Ok::<(), crate::TrajectoryError>(())
    })
    .expect("stream");
    assert_eq!(count, 1);
    assert_eq!(context.reserved_bytes(), 32768);
    drop(reader);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn cancellation_precedes_opening_the_input() {
    let context = context(32768);
    context.cancellation().cancel();
    let result = read_trajectory_in(
        Path::new("absent.trr"),
        &TrajectoryReaderOptions {
            memory_limit_bytes: 32768,
            ..TrajectoryReaderOptions::default()
        },
        &context,
    );
    assert!(matches!(
        result,
        Err(TrajectoryIoError::Reader(TrajectoryError::Cancelled))
    ));
    assert_eq!(context.reserved_bytes(), 0);
}
