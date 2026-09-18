use super::*;
use molframe_traj::{Frame, StreamingReader};

#[test]
fn surface_frames_use_one_borrowed_buffer_and_emit_only_complete_totals() {
    let positions = [[0.0; 3], [1.5, 0.0, 0.0], [5.0; 3]];
    let radii = [1.0; 3];
    let expected: f64 =
        molframe_surface::shrake_rupley(&positions, &radii, 1.4, 64, &ExecutionContext::default())
            .expect("reference")
            .iter()
            .sum();
    for workers in [1, 2, 4, 8] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("context");
        let mut reader = StreamingReader::new(
            (0..100).map(|_| Frame {
                positions: positions.to_vec(),
            }),
            3,
        );
        let mut seen = 0;
        let frames = sasa_stream(
            &mut reader,
            &radii,
            1.4,
            64,
            &context,
            4096,
            |frame, time, area| {
                assert_eq!(frame, seen);
                assert_eq!(time, None);
                assert_eq!(area.to_bits(), expected.to_bits());
                assert_eq!(context.reserved_bytes(), 4096 + 3 * 4 + 64 * 24);
                seen += 1;
                Ok(())
            },
        )
        .expect("stream");
        assert_eq!(frames, 100);
        assert_eq!(context.reserved_bytes(), 0);
    }
}

#[test]
fn invalid_radii_do_not_advance_and_a_sink_error_stops_at_that_frame() {
    let context = ExecutionContext::default();
    let mut reader = StreamingReader::new(
        (0..3).map(|_| Frame {
            positions: vec![[0.0; 3]],
        }),
        1,
    );
    assert!(
        sasa_stream(
            &mut reader,
            &[f32::NAN],
            1.4,
            32,
            &context,
            4096,
            |_, _, _| panic!("invalid radii")
        )
        .is_err()
    );
    assert!(matches!(
        sasa_stream(
            &mut reader,
            &[1.0],
            1.4,
            32,
            &context,
            4096,
            |frame, _, _| {
                assert_eq!(frame, 0);
                Err(TrajectoryError::Cancelled.into())
            }
        ),
        Err(SasaStreamError::Trajectory(TrajectoryError::Cancelled))
    ));
    assert_eq!(context.reserved_bytes(), 0);
}
