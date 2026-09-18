use super::*;

fn fixed(value: f64, decimals: u32) -> String {
    let mut buffer = Vec::new();
    push_fixed(&mut buffer, value, decimals);
    String::from_utf8_lossy(&buffer).into_owned()
}

fn integer(value: u64) -> String {
    let mut buffer = Vec::new();
    push_u64(&mut buffer, value);
    String::from_utf8_lossy(&buffer).into_owned()
}

#[test]
fn integers_match_the_standard_formatter() {
    for value in [0, 1, 9, 10, 99, 100, 12345, u64::from(u32::MAX), u64::MAX] {
        assert_eq!(integer(value), format!("{value}"));
    }
}

#[test]
fn signed_integers_match_the_standard_formatter() {
    for value in [0_i32, 1, -1, 42, -42, i32::MAX, i32::MIN] {
        let mut buffer = Vec::new();
        push_i32(&mut buffer, value);
        assert_eq!(String::from_utf8_lossy(&buffer), format!("{value}"));
    }
}

#[test]
fn three_decimal_places_match_the_standard_formatter() {
    let values = [
        0.0, 1.0, -1.0, 0.5, -0.5, 12.345, -12.345, 999.999, -999.999, 0.001, -0.001, 1234.5678,
        -1234.5678, 100.0,
    ];
    for value in values {
        assert_eq!(
            fixed(value, 3),
            format!("{value:.3}"),
            "differs for {value}"
        );
    }
}

/// A fixed-point rendering as the integer it denotes in its last place.
///
/// Comparing renderings through `f64` would reintroduce the rounding this whole
/// module exists to control, so the comparison is done in integer space.
fn scaled(rendering: &str) -> i64 {
    let digits: String = rendering.chars().filter(|c| *c != '.').collect();
    let Ok(value) = digits.parse::<i64>() else {
        panic!("rendering is not fixed point: {rendering}")
    };
    value
}

/// The gap between a rendering and the standard formatter's, in last places.
fn last_place_gap(value: f64, decimals: u32) -> i64 {
    let Ok(width) = usize::try_from(decimals) else {
        panic!("decimal count does not fit a usize")
    };
    (scaled(&fixed(value, decimals)) - scaled(&format!("{value:.width$}"))).abs()
}

#[test]
fn an_exact_decimal_midpoint_differs_by_at_most_one_last_place() {
    for (value, decimals) in [(0.0005_f64, 3_u32), (100.005, 2), (0.125, 2), (2.675, 2)] {
        let gap = last_place_gap(value, decimals);
        assert!(gap <= 1, "gap of {gap} last places for {value}");
    }
}

#[test]
fn two_decimal_places_match_the_standard_formatter() {
    for value in [0.0, 20.0, 37.5, -37.5, 99.99, -0.004, 1.0, -1.0] {
        assert_eq!(
            fixed(value, 2),
            format!("{value:.2}"),
            "differs for {value}"
        );
    }
}

#[test]
fn a_swept_coordinate_range_matches_the_standard_formatter() {
    let mut value = -500.0_f64;
    while value < 500.0 {
        assert_eq!(
            fixed(value, 3),
            format!("{value:.3}"),
            "differs for {value}"
        );
        value += 0.037;
    }
}

#[test]
fn a_value_that_is_not_a_number_is_written_as_unknown() {
    assert_eq!(fixed(f64::NAN, 3), "?");
    assert_eq!(fixed(f64::INFINITY, 3), "?");
    assert_eq!(fixed(f64::NEG_INFINITY, 3), "?");
}

#[test]
fn a_value_beyond_the_exact_integer_range_is_written_as_unknown() {
    assert_eq!(fixed(1.0e300, 3), "?");
}

#[test]
fn zero_decimals_emits_no_point() {
    assert_eq!(fixed(42.7, 0), "43");
}

#[test]
fn a_negative_value_rounding_to_zero_keeps_its_sign() {
    assert_eq!(fixed(-0.0004, 3), format!("{:.3}", -0.0004));
    assert_eq!(fixed(-0.0004, 3), "-0.000");
}
