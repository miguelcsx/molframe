use super::{f64_to_f32, usize_to_u32};
use num_traits::ToPrimitive;

/// Compares by bit pattern: these tests are about which `f32` a conversion picks,
/// and two values that differ in any bit are different outcomes.
fn same(left: f32, right: f32) -> bool {
    left.to_bits() == right.to_bits()
}

fn nearest(value: f64) -> f32 {
    let Some(narrowed) = value.to_f32() else {
        panic!("{value} is representable as f32")
    };
    narrowed
}

#[test]
fn a_representable_value_narrows_to_the_nearest_f32() {
    for value in [0.0_f64, 1.5, -2.25, 27.34, -103.7, 1e-9, -0.0] {
        let narrowed = f64_to_f32(value);
        assert!(
            same(narrowed, nearest(value)),
            "{value} narrows to {narrowed}"
        );
    }
}

#[test]
fn a_value_beyond_the_range_saturates_instead_of_becoming_an_infinity() {
    let high = f64_to_f32(f64::MAX);
    let low = f64_to_f32(f64::MIN);
    assert!(high.is_finite(), "an over-range value must stay finite");
    assert!(low.is_finite(), "an under-range value must stay finite");
    assert!(same(high, f32::MAX), "{high} must saturate to f32::MAX");
    assert!(same(low, f32::MIN), "{low} must saturate to f32::MIN");
}

#[test]
fn recorded_non_finite_values_are_carried_through() {
    assert!(f64_to_f32(f64::NAN).is_nan());
    assert!(
        same(f64_to_f32(f64::INFINITY), f32::INFINITY),
        "a recorded infinity is not a saturation"
    );
    assert!(
        same(f64_to_f32(f64::NEG_INFINITY), f32::NEG_INFINITY),
        "a recorded negative infinity is not a saturation"
    );
}

#[test]
fn an_index_wider_than_the_range_saturates() {
    assert_eq!(usize_to_u32(0), 0);
    assert_eq!(usize_to_u32(7), 7);
    assert_eq!(usize_to_u32(usize::MAX), u32::MAX);
}
