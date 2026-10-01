use super::GaussianFormFactor;
use molframe_core::Element;

fn element(symbol: &str) -> Element {
    match Element::from_symbol(symbol) {
        Some(element) => element,
        None => panic!("{symbol} is an element"),
    }
}

#[test]
fn the_forward_amplitude_is_the_atomic_number_within_the_fit() {
    // At zero angle every neutral atom scatters its electron count; the
    // four-Gaussian fit reproduces that to a few thousandths of an electron.
    for symbol in ["H", "C", "N", "O", "S", "Fe", "Zn", "Au", "U"] {
        let element = element(symbol);
        let Some(factor) = GaussianFormFactor::xray(element) else {
            panic!("{symbol} has a form factor")
        };
        let forward = factor.value(0.0);
        let electrons = f64::from(element.atomic_number());
        assert!(
            (forward - electrons).abs() < 0.05 * electrons.max(1.0),
            "{symbol}: {forward} vs {electrons}"
        );
    }
}

#[test]
fn the_amplitude_falls_with_angle_and_matches_a_tabulated_carbon_value() {
    let Some(carbon) = GaussianFormFactor::xray(element("C")) else {
        panic!("carbon has a form factor")
    };
    let mut previous = carbon.value(0.0);
    for step in 1..=20 {
        let current = carbon.value(f64::from(step) * 0.05);
        assert!(current < previous);
        previous = current;
    }
}

#[test]
fn values_match_the_reference_implementation_to_single_precision() {
    // Reference values from Gemmi 0.7.5 at (sin(theta)/lambda)^2 = 0.09, which
    // stores the coefficients in single precision.
    for (symbol, expected) in [
        ("H", 0.331_106_54),
        ("C", 2.494_160_4),
        ("O", 4.089_334),
        ("Fe", 16.748_018),
        ("Au", 59.405_533),
    ] {
        let Some(factor) = GaussianFormFactor::xray(element(symbol)) else {
            panic!("{symbol} has a form factor")
        };
        assert!(
            (factor.value(0.09) - expected).abs() < 1e-4 * expected.max(1.0),
            "{symbol}: {} vs {expected}",
            factor.value(0.09)
        );
    }
}

#[test]
fn elements_without_a_tabulated_entry_have_none() {
    assert!(GaussianFormFactor::xray(element("Fm")).is_none());
}
