//! Explicit count conversion for audit normalisation.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    let Some(converted) = value.to_f64() else {
        return f64::INFINITY;
    };
    converted
}
