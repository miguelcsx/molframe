use super::*;
use crate::{Frame, StreamingReader};

#[test]
fn long_stream_reuses_one_buffer_and_releases_its_reservation() {
    let context = ExecutionContext::default();
    let mut reader = StreamingReader::new(
        (0..10000).map(|_| Frame {
            positions: vec![[1.0, 2.0, 3.0]; 10],
        }),
        10,
    );
    let mut pointer = None;
    let mut count = 0;
    let frames = run_analysis_stream(&mut reader, &context, 4096, |frame, _| {
        if let Some(pointer) = pointer {
            assert_eq!(frame.positions.as_ptr(), pointer);
        }
        pointer = Some(frame.positions.as_ptr());
        count += 1;
        Ok::<(), crate::TrajectoryError>(())
    })
    .expect("stream");
    assert_eq!(frames, 10000);
    assert_eq!(count, 10000);
    assert_eq!(context.peak_reserved_bytes(), 4096);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn oversized_frame_is_refused_without_consuming_the_source() {
    let context = ExecutionContext::default();
    let mut reader = StreamingReader::new(
        std::iter::once(Frame {
            positions: vec![[1.0; 3]; 100],
        }),
        100,
    );
    assert!(matches!(
        run_analysis_stream(&mut reader, &context, 100, |_, _| Ok::<
            (),
            crate::TrajectoryError,
        >(())),
        Err(TrajectoryError::MemoryLimit { .. })
    ));
    let mut found = false;
    run_analysis_stream(&mut reader, &context, 4096, |frame, _| {
        found = true;
        assert_eq!(frame.frame, 0);
        Ok::<(), crate::TrajectoryError>(())
    })
    .expect("retry");
    assert!(found);
}
