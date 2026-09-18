use super::*;
use crate::{write_dlpoly_config, write_dlpoly_history};

const CONFIG: &str = "periodic state\n 2 3 2\n\
10 0 0\n1 9 0\n0.5 0.25 8\n\
C 1\n1 2 3\n0.1 0.2 0.3\n100 200 300\n\
O 2\n4 5 6\n0.4 0.5 0.6\n-100 -200 -300\n";

const HISTORY: &str = "trajectory\n2 3 2\n\
timestep 10 2 2 3 0.5\n10 0 0\n1 9 0\n0.5 0.25 8\n\
C 1 12.011 0.2\n1 2 3\n0.1 0.2 0.3\n100 200 300\n\
O 2 15.999 -0.2\n4 5 6\n0.4 0.5 0.6\n-100 -200 -300\n\
timestep 20 2 2 3 1.0\n10 0 0\n1 9 0\n0.5 0.25 8\n\
C 1 12.011 0.2\n2 3 4\n0.2 0.3 0.4\n200 300 400\n\
O 2 15.999 -0.2\n5 6 7\n0.5 0.6 0.7\n-200 -300 -400\n";

#[test]
fn config_round_trips_cell_velocities_and_canonical_forces() {
    let config =
        parse_dlpoly_config(CONFIG).unwrap_or_else(|error| panic!("CONFIG parse failed: {error}"));
    assert_eq!(
        config.frame.forces.as_ref().map(|forces| forces[0]),
        Some([1.0, 2.0, 3.0])
    );
    let text =
        write_dlpoly_config(&config).unwrap_or_else(|error| panic!("CONFIG write failed: {error}"));
    let observed =
        parse_dlpoly_config(&text).unwrap_or_else(|error| panic!("written CONFIG failed: {error}"));
    assert_eq!(observed, config);
}

#[test]
fn history_round_trips_frames_and_reuses_one_topology() {
    let history = parse_dlpoly_history(HISTORY)
        .unwrap_or_else(|error| panic!("HISTORY parse failed: {error}"));
    assert_eq!(history.frames.len(), 2);
    assert_eq!(history.atoms.len(), 2);
    assert_eq!(history.frames[1].to_timestep(1).time, Some(1.0));
    let text = write_dlpoly_history(&history)
        .unwrap_or_else(|error| panic!("HISTORY write failed: {error}"));
    let observed = parse_dlpoly_history(&text)
        .unwrap_or_else(|error| panic!("written HISTORY failed: {error}"));
    assert_eq!(observed, history);
}

#[test]
fn history_rejects_topology_drift() {
    let changed = HISTORY.replacen("O 2 15.999 -0.2", "N 2 15.999 -0.2", 1);
    let changed = changed.replacen("O 2 15.999 -0.2", "S 2 15.999 -0.2", 1);
    assert_eq!(
        parse_dlpoly_history(&changed),
        Err(DlPolyError::TopologyDrift)
    );
}
