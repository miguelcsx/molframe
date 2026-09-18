use super::PyTrajectory;
use molframe::traj::{FormatMetadata, Timestep, TrajectoryFormat, TrajectoryMetadata};
use pyo3::Python;

#[test]
fn normalized_storage_preserves_native_streams_and_metadata() {
    let frames = vec![
        Timestep {
            frame: 0,
            time: Some(1.0),
            positions: vec![[1.0, 2.0, 3.0]],
            velocities: Some(vec![[0.1, 0.2, 0.3]]),
            ..Timestep::default()
        },
        Timestep {
            frame: 1,
            time: Some(2.0),
            positions: vec![[4.0, 5.0, 6.0]],
            velocities: Some(vec![[0.4, 0.5, 0.6]]),
            ..Timestep::default()
        },
    ];
    let metadata = TrajectoryMetadata {
        steps: Some(vec![10, 20]),
        format: FormatMetadata::Xtc {
            precision: vec![1_000.0, 1_000.0],
        },
    };
    let trajectory = PyTrajectory::from_frames(&frames, Some(TrajectoryFormat::Xtc), metadata)
        .expect("normalization should succeed");
    let restored = Python::attach(|py| trajectory.to_data(py, TrajectoryFormat::Xtc))
        .expect("materialization should succeed");
    assert_eq!(restored.frames, frames);
    assert_eq!(restored.metadata.steps, Some(vec![10, 20]));
}

#[test]
fn mixed_optional_streams_are_rejected_without_sentinels() {
    let frames = vec![
        Timestep {
            positions: vec![[1.0, 2.0, 3.0]],
            time: Some(1.0),
            ..Timestep::default()
        },
        Timestep {
            positions: vec![[4.0, 5.0, 6.0]],
            ..Timestep::default()
        },
    ];
    assert!(PyTrajectory::from_frames(&frames, None, TrajectoryMetadata::default()).is_err());
}
