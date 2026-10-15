use super::sugar_pucker;

#[test]
fn equal_torsions_have_no_pucker() {
    let pucker = sugar_pucker([5.0, 5.0, 5.0, 5.0, 5.0]);
    assert!(pucker.phase_degrees.abs() < 1e-9);
    assert!(pucker.amplitude.abs() < 1e-9);
}

#[test]
fn sweep_recovers_phase_and_amplitude_including_90_and_270() {
    for amplitude in [10.0, 38.0] {
        for p in 0..360 {
            let phase = f64::from(p);
            let mut nu = [0.0; 5];
            for (j, slot) in nu.iter_mut().enumerate() {
                let angle = phase + 144.0 * (f64::from(u8::try_from(j).unwrap_or(0)) - 2.0);
                *slot = amplitude * angle.to_radians().cos();
            }
            let pucker = sugar_pucker(nu);
            let mut diff = (pucker.phase_degrees - phase).abs();
            diff = diff.min(360.0 - diff);
            assert!(diff < 1e-6, "P {phase}: got {}", pucker.phase_degrees);
            assert!(
                (pucker.amplitude - amplitude).abs() < 1e-9,
                "P {phase}: amplitude {}",
                pucker.amplitude
            );
        }
    }
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
