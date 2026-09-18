//! Counting past the end of a fixed-width field.
//!
//! Five columns hold serial numbers up to 99999 and four hold residue numbers up
//! to 9999, which large assemblies passed long ago. Hybrid-36 continues the
//! sequence in base 36 rather than overflowing: once the decimal range is
//! exhausted the field switches to upper-case alphanumerics, then to lower-case.
//!
//! Reading it is unconditional, because a file that uses it is otherwise
//! unreadable. Writing it is not, because a consumer that does not implement it
//! would silently read the wrong number — so it is offered and never assumed.

/// The digits of each case, in value order.
const UPPER: &[u8; 36] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const LOWER: &[u8; 36] = b"0123456789abcdefghijklmnopqrstuvwxyz";

/// Decodes a field that may be decimal or hybrid-36.
///
/// # Examples
///
/// ```
/// use molframe_pdb::hybrid36::decode;
///
/// assert_eq!(decode("99999", 5), Some(99_999));
/// assert_eq!(decode("A0000", 5), Some(100_000));
/// assert_eq!(decode("  -12", 5), Some(-12));
/// ```
#[must_use]
pub fn decode(field: &str, width: u32) -> Option<i64> {
    let trimmed = field.trim();
    if trimmed.is_empty() {
        return None;
    }
    if trimmed
        .bytes()
        .all(|byte| byte.is_ascii_digit() || byte == b'-' || byte == b'+')
    {
        return trimmed.parse().ok();
    }
    if trimmed.len() != width as usize {
        return None;
    }

    let decimal_range = 10i64.checked_pow(width)?;
    let case_span = case_span(width)?;
    let lower = trimmed.bytes().any(|byte| byte.is_ascii_lowercase());
    let value = from_base36(trimmed, lower)?;
    let start = case_start(width)?;

    let offset = value.checked_sub(start)?;
    if offset < 0 || offset >= case_span {
        return None;
    }
    if lower {
        Some(decimal_range + case_span + offset)
    } else {
        Some(decimal_range + offset)
    }
}

/// Encodes a value into `width` columns, using hybrid-36 only when it must.
///
/// Returns nothing when the value does not fit even in hybrid-36, because a
/// number that cannot be written is not a number to write approximately.
#[must_use]
pub fn encode(value: i64, width: u32) -> Option<String> {
    let columns = width as usize;
    let decimal_range = 10i64.checked_pow(width)?;

    if value < 0 {
        let text = value.to_string();
        return (text.len() <= columns).then(|| format!("{text:>columns$}"));
    }
    if value < decimal_range {
        return Some(format!("{value:>columns$}"));
    }

    let case_span = case_span(width)?;
    let offset = value - decimal_range;
    let lower = offset >= case_span;
    let offset = if lower { offset - case_span } else { offset };
    if offset >= case_span {
        return None;
    }
    to_base36(case_start(width)? + offset, width, lower)
}

/// Returns true when a value needs hybrid-36 to fit `width` columns.
#[must_use]
pub fn needs_encoding(value: i64, width: u32) -> bool {
    10i64.checked_pow(width).is_some_and(|range| value >= range)
}

/// The first base-36 value of a case range: `A000…` read in that case's digits.
///
/// Both cases start at the digit worth ten followed by zeroes, because each case
/// has its own alphabet; which alphabet was used is what tells the two ranges
/// apart, not the value.
fn case_start(width: u32) -> Option<i64> {
    10i64.checked_mul(36i64.checked_pow(width.checked_sub(1)?)?)
}

/// How many values each case range holds.
fn case_span(width: u32) -> Option<i64> {
    36i64.checked_pow(width)?.checked_sub(case_start(width)?)
}

fn from_base36(text: &str, lower: bool) -> Option<i64> {
    let alphabet = if lower { LOWER } else { UPPER };
    let mut value = 0i64;
    for byte in text.bytes() {
        let digit = alphabet.iter().position(|candidate| *candidate == byte)?;
        value = value
            .checked_mul(36)?
            .checked_add(i64::try_from(digit).ok()?)?;
    }
    Some(value)
}

fn to_base36(mut value: i64, width: u32, lower: bool) -> Option<String> {
    let alphabet = if lower { LOWER } else { UPPER };
    let mut digits = vec![b'0'; usize::try_from(width).ok()?];
    for slot in digits.iter_mut().rev() {
        *slot = *alphabet.get(usize::try_from(value % 36).ok()?)?;
        value /= 36;
    }
    if value != 0 {
        return None;
    }
    String::from_utf8(digits).ok()
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
