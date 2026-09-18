use super::{Frame, Trajectory};

fn frame(x: f32) -> Frame {
    Frame {
        positions: vec![[x, 0.0, 0.0]],
    }
}

#[test]
fn frames_are_reachable_by_index() {
    let trajectory = Trajectory::from_frames(vec![frame(0.0), frame(1.0), frame(2.0)])
        .expect("fixed-width trajectory");
    assert_eq!(trajectory.len(), 3);
    assert_eq!(
        trajectory.frame(1).map(|frame| frame.positions),
        Some(&[[1.0, 0.0, 0.0]][..])
    );
    assert!(trajectory.frame(5).is_none());
}

#[test]
fn slicing_selects_and_reorders_frames() {
    let trajectory = Trajectory::from_frames(vec![frame(0.0), frame(1.0), frame(2.0)])
        .expect("fixed-width trajectory");
    let sliced = trajectory.slice(&[2, 0]).expect("valid slice dimensions");
    assert_eq!(sliced.len(), 2);
    assert_eq!(
        sliced.frame(0).map(|frame| frame.positions),
        Some(&[[2.0, 0.0, 0.0]][..])
    );
    assert_eq!(
        sliced.frame(1).map(|frame| frame.positions),
        Some(&[[0.0, 0.0, 0.0]][..])
    );
}

#[test]
fn slicing_skips_out_of_range_indices() {
    let trajectory = Trajectory::from_frames(vec![frame(0.0)]).expect("fixed-width trajectory");
    assert_eq!(
        trajectory
            .slice(&[0, 9])
            .expect("valid slice dimensions")
            .len(),
        1
    );
}

#[test]
fn an_empty_trajectory_reports_empty() {
    assert!(Trajectory::default().is_empty());
}
