use super::{TrajectoryInterpolation, TrajectoryInterpolationError, interpolate_trajectory_frames};

#[test]
fn linear_interpolation_hits_the_central_endpoints() {
    let f0 = [[-1.0, 0.0, 0.0]];
    let f1 = [[0.0, 0.0, 0.0]];
    let f2 = [[2.0, 0.0, 0.0]];
    let f3 = [[3.0, 0.0, 0.0]];
    let frames = [&f0[..], &f1[..], &f2[..], &f3[..]];
    assert_eq!(
        interpolate_trajectory_frames(frames, 0.0, TrajectoryInterpolation::Linear),
        Ok(vec![f1[0]])
    );
    assert_eq!(
        interpolate_trajectory_frames(frames, 1.0, TrajectoryInterpolation::Linear),
        Ok(vec![f2[0]])
    );
}

#[test]
fn centripetal_interpolation_is_finite_and_deterministic() {
    let f0 = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let f1 = [[1.0, 0.0, 0.0], [1.0, 1.5, 0.0]];
    let f2 = [[2.0, 1.0, 0.0], [2.0, 2.0, 0.5]];
    let f3 = [[4.0, 1.0, 0.0], [3.0, 3.0, 1.0]];
    let frames = [&f0[..], &f1[..], &f2[..], &f3[..]];
    let left =
        interpolate_trajectory_frames(frames, 0.5, TrajectoryInterpolation::CentripetalCatmullRom);
    let right =
        interpolate_trajectory_frames(frames, 0.5, TrajectoryInterpolation::CentripetalCatmullRom);
    assert_eq!(left, right);
    let Ok(points) = left else {
        panic!("cubic interpolation should succeed")
    };
    assert!(points.iter().flatten().all(|value| value.is_finite()));
}

#[test]
fn mismatched_topology_is_rejected() {
    let one = [[0.0; 3]];
    let two = [[0.0; 3], [1.0; 3]];
    let result = interpolate_trajectory_frames(
        [&one, &one, &two, &one],
        0.5,
        TrajectoryInterpolation::Linear,
    );
    assert_eq!(result, Err(TrajectoryInterpolationError::TopologyMismatch));
}

#[test]
fn repeated_cubic_knots_fall_back_to_linear_motion() {
    let repeated = [[1.0, 0.0, 0.0]];
    let end = [[3.0, 0.0, 0.0]];
    let result = interpolate_trajectory_frames(
        [&repeated, &repeated, &end, &end],
        0.25,
        TrajectoryInterpolation::CentripetalCatmullRom,
    );
    assert_eq!(result, Ok(vec![[1.5, 0.0, 0.0]]));
}
