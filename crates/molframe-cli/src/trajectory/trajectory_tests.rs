use super::command::apply_stride;
use super::summary::Summary;
use molframe::traj::{
    FormatMetadata, Timestep, TrajectoryData, TrajectoryFormat, TrajectoryMetadata,
};

#[test]
fn summary_reports_optional_stream_coverage() {
    let frames = vec![
        Timestep {
            time: Some(0.0),
            positions: vec![[1.0, 2.0, 3.0]],
            velocities: Some(vec![[0.0; 3]]),
            ..Timestep::default()
        },
        Timestep {
            time: Some(2.0),
            positions: vec![[4.0, 5.0, 6.0]],
            forces: Some(vec![[0.0; 3]]),
            ..Timestep::default()
        },
    ];
    let data = TrajectoryData {
        format: TrajectoryFormat::Trr,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::None,
        },
    };
    let summary = Summary::from_data(&data);
    assert_eq!(summary.format, "trr");
    assert_eq!(summary.frames, 2);
    assert_eq!(summary.atoms, 1);
    assert_eq!(summary.first_time, Some(0.0));
    assert_eq!(summary.last_time, Some(2.0));
    assert_eq!(summary.frames_with_velocities, 1);
    assert_eq!(summary.frames_with_forces, 1);
}

#[test]
fn stride_keeps_frames_and_steps_aligned() {
    let mut data = TrajectoryData {
        format: TrajectoryFormat::Xtc,
        frames: [0.0_f32, 1.0, 2.0, 3.0, 4.0]
            .into_iter()
            .enumerate()
            .map(|(frame, coordinate)| Timestep {
                frame,
                positions: vec![[coordinate, 0.0, 0.0]],
                ..Timestep::default()
            })
            .collect(),
        metadata: TrajectoryMetadata {
            steps: Some(vec![10, 20, 30, 40, 50]),
            format: FormatMetadata::None,
        },
    };
    apply_stride(&mut data, 2);
    assert_eq!(data.frames.len(), 3);
    assert_eq!(data.frames[1].frame, 2);
    assert_eq!(data.metadata.steps, Some(vec![10, 30, 50]));
}
