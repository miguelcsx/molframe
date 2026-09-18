//! Checked numeric conversions shared by geometric kernels.

use num_traits::ToPrimitive;

/// Largest consecutive integer representable exactly by an IEEE-754 `f64`.
const MAX_EXACT_INTEGER_F64: u64 = 1_u64 << f64::MANTISSA_DIGITS;

pub(crate) fn exact_count(value: usize) -> Option<f64> {
    exact_u64(u64::try_from(value).ok()?)
}

pub(crate) fn exact_u64(value: u64) -> Option<f64> {
    if value > MAX_EXACT_INTEGER_F64 {
        return None;
    }
    value.to_f64()
}
