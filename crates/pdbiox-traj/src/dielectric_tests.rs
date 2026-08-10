use super::{DielectricOptions, dielectric_from_dipoles};

#[test]
fn dipole_fluctuation_and_prefactor_are_reported_separately() {
    let dipoles = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let Ok(estimate) = dielectric_from_dipoles(
        &dipoles,
        DielectricOptions {
            fluctuation_prefactor: 2.0,
        },
    ) else {
        panic!("valid dipoles");
    };
    estimate
        .mean_dipole
        .into_iter()
        .for_each(|value| assert!(value.abs() < 1.0e-12));
    assert!((estimate.fluctuation - 1.0).abs() < 1.0e-12);
    assert!((estimate.relative_permittivity - 3.0).abs() < 1.0e-12);
}

#[test]
fn no_temperature_volume_or_unit_fallback_is_used() {
    assert!(
        dielectric_from_dipoles(
            &[[0.0; 3]],
            DielectricOptions {
                fluctuation_prefactor: f64::NAN,
            },
        )
        .is_err()
    );
}
