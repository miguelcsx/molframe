//! Explicit count conversion for floating-point score normalisation.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    finite_or_infinity(value.to_f64())
}

pub(crate) fn u64_to_f64(value: u64) -> f64 {
    finite_or_infinity(value.to_f64())
}

fn finite_or_infinity(value: Option<f64>) -> f64 {
    let Some(value) = value else {
        return f64::INFINITY;
    };
    value
}
