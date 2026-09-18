use super::sugar_pucker;

#[test]
fn equal_torsions_have_zero_phase_and_amplitude_equal_to_nu2() {
    let pucker = sugar_pucker([5.0, 5.0, 5.0, 5.0, 5.0]);
    assert!(
        pucker.phase_degrees.abs() < 1e-9,
        "phase {}",
        pucker.phase_degrees
    );
    assert!((pucker.amplitude - 5.0).abs() < 1e-9);
}

#[test]
fn a_positive_out_of_plane_pattern_gives_a_quarter_turn_phase() {
    // Numerator (ν4+ν1)−(ν3+ν0) = 2, denominator 2·ν2·… = 0, so atan2 gives 90°.
    let pucker = sugar_pucker([0.0, 1.0, 0.0, 0.0, 1.0]);
    assert!(
        (pucker.phase_degrees - 90.0).abs() < 1e-9,
        "phase {}",
        pucker.phase_degrees
    );
}

#[test]
fn a_negative_numerator_wraps_into_the_full_circle() {
    // Numerator = −1, denominator 0 → atan2 = −90°, wrapped to 270°.
    let pucker = sugar_pucker([0.0, 0.0, 0.0, 1.0, 0.0]);
    assert!(
        (pucker.phase_degrees - 270.0).abs() < 1e-9,
        "phase {}",
        pucker.phase_degrees
    );
}
