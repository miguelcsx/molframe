use super::MinimalTopology;
use crate::{Frame, MemoryReader, Trajectory, TrajectoryError};

#[test]
fn atom_count_is_checked_before_a_coordinate_frame_is_consumed() {
    let trajectory = Trajectory::from_frames(vec![Frame {
        positions: vec![[0.0; 3]; 3],
    }]);
    let reader = MemoryReader::new(&trajectory);
    assert_eq!(MinimalTopology::new(3).validate_reader(&reader), Ok(()));
    assert_eq!(
        MinimalTopology::new(4).validate_reader(&reader),
        Err(TrajectoryError::AtomCountMismatch {
            expected: 4,
            found: 3,
        })
    );
}
