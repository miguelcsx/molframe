//! Explicit saturating conversions for compact chemistry values.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    let Ok(converted) = u32::try_from(value) else {
        return u32::MAX;
    };
    converted
}

/// Narrows a coordinate to `f32`, saturating at the ends of the range.
///
/// A structure stores coordinates as `f32`, so a site whose value cannot be
/// represented is a property of the input rather than a rounding decision this
/// layer gets to make. A cast would turn an out-of-range magnitude into an
/// infinity silently; here it saturates into the finite range instead.
///
/// NaN and the infinities are carried through as themselves — "not a number" and
/// "not finite" are things a reader may have recorded deliberately, and they are
/// distinguishable from a finite value that merely exceeded the range.
pub(crate) fn f64_to_f32(value: f64) -> f32 {
    if value.is_nan() {
        return f32::NAN;
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            f32::NEG_INFINITY
        } else {
            f32::INFINITY
        };
    }
    let Some(narrowed) = value.to_f32() else {
        return saturate(value);
    };
    // A finite `f64` outside the `f32` range narrows to an infinity, which is the
    // point the conversion failed rather than a value the input held.
    if narrowed.is_infinite() {
        return saturate(value);
    }
    narrowed
}

/// The finite end of the `f32` range nearest to a finite `f64` outside it.
fn saturate(value: f64) -> f32 {
    if value.is_sign_negative() {
        f32::MIN
    } else {
        f32::MAX
    }
}

#[cfg(test)]
#[path = "numeric_tests.rs"]
mod tests;
