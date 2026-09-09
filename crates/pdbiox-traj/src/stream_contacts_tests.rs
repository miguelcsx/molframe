use super::*;
use crate::{Frame, StreamingReader};

#[test]
fn long_contact_stream_releases_each_index_and_preserves_frame_order() {
    let context = ExecutionContext::builder()
        .worker_budget(4)
        .build()
        .expect("context");
    let mut reader = StreamingReader::new(
        (0..1000).map(|_| Frame {
            positions: vec![[0.0; 3], [1.0, 0.0, 0.0], [5.0; 3]],
        }),
        3,
    );
    let mut seen = 0;
    let frames = contact_counts_stream(
        &mut reader,
        1.0,
        SpatialSearchOptions::with_backend(pdbiox_spatial::SpatialBackend::CellList),
        &context,
        4096,
        |index, time, count| {
            assert_eq!(index, seen);
            assert_eq!(time, None);
            assert_eq!(count, 1);
            // The index and partial reductions have already been released.
            assert_eq!(context.reserved_bytes(), 4096);
            seen += 1;
            Ok(())
        },
    )
    .expect("stream");
    assert_eq!(frames, 1000);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn invalid_cutoff_does_not_consume_a_frame_and_sink_errors_stop_the_stream() {
    let context = ExecutionContext::default();
    let mut reader = StreamingReader::new(
        (0..3).map(|_| Frame {
            positions: vec![[0.0; 3]; 2],
        }),
        2,
    );
    assert!(
        contact_counts_stream(
            &mut reader,
            -1.0,
            SpatialSearchOptions::default(),
            &context,
            4096,
            |_, _, _| panic!("invalid input emitted")
        )
        .is_err()
    );
    let error = contact_counts_stream(
        &mut reader,
        1.0,
        SpatialSearchOptions::default(),
        &context,
        4096,
        |index, _, count| {
            assert_eq!(index, 0);
            assert_eq!(count, 1);
            Err(TrajectoryError::Cancelled)
        },
    )
    .expect_err("sink failure");
    assert_eq!(error, TrajectoryError::Cancelled);
    assert_eq!(context.reserved_bytes(), 0);
}
