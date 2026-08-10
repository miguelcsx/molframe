//! Explicit count conversion for floating-point score normalisation.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::INFINITY,
    }
}

pub(crate) fn u64_to_f64(value: u64) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::INFINITY,
    }
}
