use super::{SurvivalMode, water_dynamics};

#[test]
fn continuous_and_intermittent_survival_differ_after_a_gap() {
    let occupancy = vec![vec![true], vec![false], vec![true]];
    let Ok(continuous) = water_dynamics(&occupancy, 2, 0.5, SurvivalMode::Continuous) else {
        panic!("valid occupancy");
    };
    let Ok(intermittent) = water_dynamics(&occupancy, 2, 0.5, SurvivalMode::Intermittent) else {
        panic!("valid occupancy");
    };
    assert_eq!(continuous.survival[2].probability, Some(0.0));
    assert_eq!(intermittent.survival[2].probability, Some(1.0));
}

#[test]
fn residence_time_integrates_the_explicit_frame_duration() {
    let occupancy = vec![vec![true], vec![true], vec![true]];
    let Ok(result) = water_dynamics(&occupancy, 2, 2.0, SurvivalMode::Continuous) else {
        panic!("valid occupancy");
    };
    assert_eq!(result.residence_time, Some(4.0));
}

#[test]
fn missing_origins_are_reported_as_absence_not_zero_probability() {
    let occupancy = vec![vec![false], vec![false]];
    let Ok(result) = water_dynamics(&occupancy, 1, 1.0, SurvivalMode::Intermittent) else {
        panic!("valid empty occupancy state");
    };
    assert_eq!(result.survival[1].probability, None);
    assert_eq!(result.residence_time, None);
}
