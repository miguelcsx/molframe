use super::*;
use crate::{Frame, StreamingReader};

fn context(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("context")
}

fn reader(frames: usize) -> StreamingReader<impl Iterator<Item = Frame>> {
    StreamingReader::new(
        (0..frames).map(|i| Frame {
            positions: vec![[f32::from(u16::try_from(i).expect("small")), 0.0, 0.0]; 3],
        }),
        3,
    )
}

fn collect(workers: usize, frames: usize) -> Vec<(usize, f32)> {
    let context = context(workers);
    let mut out = Vec::new();
    let count = pipeline::<_, _, TrajectoryError, _, _>(
        &mut reader(frames),
        &context,
        4096,
        Duration::ZERO,
        |frame, _| Ok(frame.positions[0][0] * 2.0),
        |frame, value| {
            out.push((frame.frame, value));
            Ok(())
        },
    )
    .expect("pipeline");
    assert_eq!(count, frames as u64);
    assert_eq!(context.reserved_bytes(), 0);
    out
}

#[test]
fn results_are_ordered_and_identical_at_every_worker_count() {
    for frames in [0, 1, 5, 64, 65, 300] {
        let serial = collect(1, frames);
        assert_eq!(serial.len(), frames);
        for workers in [2, 3, 8] {
            assert_eq!(collect(workers, frames), serial, "{workers}/{frames}");
        }
    }
}

#[test]
fn compute_failure_emits_earlier_frames_then_stops() {
    let context = context(4);
    let mut seen = Vec::new();
    let error = pipeline(
        &mut reader(100),
        &context,
        4096,
        Duration::ZERO,
        |frame, _| {
            if frame.frame == 7 {
                Err(TrajectoryError::DegenerateFit)
            } else {
                Ok(frame.frame)
            }
        },
        |_, value| {
            seen.push(value);
            Ok(())
        },
    )
    .expect_err("failure");
    assert_eq!(error, TrajectoryError::DegenerateFit);
    assert_eq!(seen, (0..7).collect::<Vec<_>>());
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn reservation_stays_within_two_windows() {
    let context = context(2);
    let mut peak = 0;
    pipeline::<_, _, TrajectoryError, _, _>(
        &mut reader(200),
        &context,
        4096,
        Duration::ZERO,
        |_, _| Ok(()),
        |_, ()| {
            peak = peak.max(context.reserved_bytes());
            Ok(())
        },
    )
    .expect("pipeline");
    assert_eq!(peak, 2 * 8 * 4096);
}

#[test]
fn cancellation_stops_reading() {
    let context = context(4);
    context.cancellation().cancel();
    let error = run_frame_pipeline::<_, _, TrajectoryError, _, _>(
        &mut reader(10),
        &context,
        4096,
        |_, _| Ok(()),
        |_, ()| panic!("emitted after cancel"),
    )
    .expect_err("cancelled");
    assert_eq!(error, TrajectoryError::Cancelled);
}

#[test]
fn cheap_frames_stay_serial_and_match_the_windowed_result() {
    let context = context(4);
    let mut peak = 0;
    let mut values = Vec::new();
    run_frame_pipeline::<_, _, TrajectoryError, _, _>(
        &mut reader(300),
        &context,
        4096,
        |frame, _| Ok(frame.frame),
        |_, value| {
            peak = peak.max(context.reserved_bytes());
            values.push(value);
            Ok(())
        },
    )
    .expect("serial");
    assert_eq!(values, (0..300).collect::<Vec<_>>());
    assert_eq!(peak, 4096);
}
