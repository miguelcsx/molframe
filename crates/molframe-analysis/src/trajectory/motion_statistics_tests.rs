use super::*;

fn ids(values: &[u32]) -> AtomSelection {
    AtomSelection::from_sorted(values.to_vec())
}

#[test]
fn intermittent_survival_and_continuous_residence_are_distinct() {
    let frames = [ids(&[1]), ids(&[]), ids(&[1])];
    let Ok(result) = water_dynamics(&frames, WaterDynamicsOptions { maximum_lag: 2 }) else {
        panic!("valid stable-ID occupancy");
    };
    assert_eq!(result[2].survival_probability, Some(1.0));
    assert_eq!(result[2].residence_probability, Some(0.0));
}

#[test]
fn dielectric_uses_the_explicit_unit_prefactor() {
    let Ok(result) = dielectric_from_dipoles(
        &[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        DielectricOptions {
            volume: 2.0,
            temperature: 2.0,
            fluctuation_prefactor: 4.0,
        },
    ) else {
        panic!("valid dipole series");
    };
    assert_eq!(result.fluctuation.to_bits(), 1.0_f64.to_bits());
    assert_eq!(result.relative_permittivity.to_bits(), 2.0_f64.to_bits());
}
