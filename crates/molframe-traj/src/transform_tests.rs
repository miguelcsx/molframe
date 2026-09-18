use super::{Center, Fit, PipelineReader, RigidTransform};
use crate::{Frame, StreamingReader, Timestep, TrajectoryReader};

#[test]
fn transforms_compose_in_declaration_order_during_streaming() {
    let frames = [Frame {
        positions: vec![[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
    }];
    let source = StreamingReader::new(frames.into_iter(), 2);
    let mut pipeline = PipelineReader::new(source)
        .then(RigidTransform::translation([1.0, 0.0, 0.0]))
        .then(Center::geometric([0, 1], [0.0; 3]));
    let mut timestep = Timestep::default();
    assert!(pipeline.read_next(&mut timestep).is_ok_and(|read| read));
    assert!((timestep.positions[0][0] + 1.0).abs() < 1.0e-6);
    assert!((timestep.positions[1][0] - 1.0).abs() < 1.0e-6);
    assert_eq!(
        pipeline.transform_names().collect::<Vec<_>>(),
        ["rigid", "center_geometric"]
    );
}

#[test]
fn fit_uses_the_shared_rigid_superposition() {
    let frames = [Frame {
        positions: vec![[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [5.0, 1.0, 0.0]],
    }];
    let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let source = StreamingReader::new(frames.into_iter(), 3);
    let mut pipeline = PipelineReader::new(source).then(Fit::new([0, 1, 2], reference));
    let mut timestep = Timestep::default();
    assert!(pipeline.read_next(&mut timestep).is_ok_and(|read| read));
    assert!(
        timestep
            .positions
            .iter()
            .zip(reference)
            .all(|(observed, expected)| observed
                .iter()
                .zip(expected)
                .all(|(left, right)| (left - right).abs() < 1.0e-5))
    );
}
