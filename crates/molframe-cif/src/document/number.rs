//! Reading the fixed-point decimals a coordinate loop is mostly made of.

use num_traits::ToPrimitive;

/// Reads `[-]digits.digits` by integer arithmetic, or `None` for anything else.
///
/// Coordinates, occupancies and temperature factors are fixed-point decimals
/// with a handful of digits, and they are most of the numbers in a coordinate
/// loop. The general float parser handles exponents, infinities and
/// seventeen-digit mantissas, and pays for that generality on every value.
///
/// The result is the same. Both the mantissa and the power of ten are exactly
/// representable at these magnitudes, so the quotient is correctly rounded, and
/// a correctly rounded quotient of the exact decimal is the same `f64` the
/// general parser produces. The bounds below are what keep that true, and the
/// tests hold the two against each other over a swept range.
pub(super) fn fixed_point(text: &str) -> Option<f64> {
    /// Largest digit count whose mantissa stays exactly representable.
    const MANTISSA_DIGITS: u32 = 15;
    /// Largest power of ten that is exactly representable.
    const SCALE_DIGITS: u32 = 22;

    let bytes = text.as_bytes();
    let (negative, digits) = match bytes.split_first() {
        Some((b'-', rest)) => (true, rest),
        Some((b'+', rest)) => (false, rest),
        _ => (false, bytes),
    };
    if digits.is_empty() {
        return None;
    }

    let mut mantissa = 0_u64;
    let mut count = 0_u32;
    let mut decimals: Option<u32> = None;

    for &byte in digits {
        match byte {
            b'0'..=b'9' => {
                count += 1;
                if count > MANTISSA_DIGITS {
                    return None;
                }
                mantissa = mantissa * 10 + u64::from(byte - b'0');
                if let Some(places) = decimals.as_mut() {
                    *places += 1;
                }
            }
            b'.' if decimals.is_none() => decimals = Some(0),
            _ => return None,
        }
    }

    let places = decimals?;
    if count == 0 || places > SCALE_DIGITS {
        return None;
    }

    let scale = 10_u64.checked_pow(places)?;
    let value = exact(mantissa)? / exact(scale)?;
    Some(if negative { -value } else { value })
}

/// A non-negative integer as an exactly representable `f64`.
///
/// Refuses anything past the mantissa's exact range, because the whole
/// argument for this path is that no rounding happens before the division.
fn exact(value: u64) -> Option<f64> {
    const EXACT_LIMIT: u64 = 1_u64 << 53;
    if value >= EXACT_LIMIT {
        return None;
    }
    value.to_f64()
}

#[cfg(test)]
#[path = "number_tests.rs"]
mod tests;
