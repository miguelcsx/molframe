//! IEEE binary16 conversion with round-to-nearest-even, kept local so the
//! dense confidence storage needs no external half-float dependency.

/// Converts a single-precision value to binary16 bits.
///
/// Values beyond the binary16 range become signed infinity; NaN stays NaN.
pub(crate) const fn f32_to_f16_bits(value: f32) -> u16 {
    let bits = value.to_bits();
    let sign = ((bits >> 16) & 0x8000) as u16;
    let exponent = ((bits >> 23) & 0xff) as i32;
    let mantissa = bits & 0x007f_ffff;
    if exponent == 0xff {
        let payload = if mantissa == 0 {
            0
        } else {
            0x0200 | (mantissa >> 13) as u16
        };
        return sign | 0x7c00 | payload;
    }
    let rebiased = exponent - 127 + 15;
    if rebiased >= 31 {
        return sign | 0x7c00;
    }
    if rebiased <= 0 {
        if rebiased < -10 {
            return sign;
        }
        let significand = mantissa | 0x0080_0000;
        let shift = (14 - rebiased) as u32;
        let kept = significand >> shift;
        let remainder = significand & ((1 << shift) - 1);
        let halfway = 1 << (shift - 1);
        let rounded = if remainder > halfway || (remainder == halfway && kept & 1 == 1) {
            kept + 1
        } else {
            kept
        };
        return sign | rounded as u16;
    }
    let combined = ((rebiased as u32) << 10) | (mantissa >> 13);
    let remainder = mantissa & 0x1fff;
    let rounded = if remainder > 0x1000 || (remainder == 0x1000 && combined & 1 == 1) {
        combined + 1
    } else {
        combined
    };
    sign | rounded as u16
}

/// Converts binary16 bits to a single-precision value exactly.
pub(crate) fn f16_bits_to_f32(bits: u16) -> f32 {
    let sign = u32::from(bits & 0x8000) << 16;
    let exponent = u32::from((bits >> 10) & 0x1f);
    let mantissa = u32::from(bits & 0x03ff);
    let pattern = match exponent {
        0 if mantissa == 0 => sign,
        0 => {
            let magnitude = f32::from(bits & 0x03ff) * (1.0 / 16_777_216.0);
            magnitude.to_bits() | sign
        }
        31 => sign | 0x7f80_0000 | (mantissa << 13),
        _ => sign | ((exponent + 112) << 23) | (mantissa << 13),
    };
    f32::from_bits(pattern)
}

#[cfg(test)]
#[path = "half_tests.rs"]
mod tests;
