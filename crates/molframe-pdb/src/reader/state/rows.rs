//! Rows that look like coordinates but do not read as one, and the two field
//! decoders every other row helper shares.
//!
//! A fixed-column file names its record in columns 1 to 6. A coordinate record
//! whose serial outgrows its five columns pushes a digit into column 6, so
//! `ATOM 100000` reads `ATOM 1` and stops naming a record this reader knows, even
//! though the keyword is plainly there. Telling that apart from a line the reader
//! simply has no interest in is the difference between a caller fixing their file
//! and a caller being told the file holds no coordinates at all.

use crate::fixed;
use crate::hybrid36;
use crate::reader::lines::Line;
use molframe_core::column::Presence;
use molframe_core::optional::OptionalI32;

/// Whether a line carries a coordinate keyword but not as its record name.
///
/// A well-formed `ATOM` or `HETATM` row has the keyword at columns 1 to 4 or 1 to
/// 6 followed by spaces. When the serial outgrows its field the digits shift left
/// into the name, so what remains after the keyword is digits rather than padding.
#[must_use]
pub(crate) fn starts_with_coordinate_keyword(line: &str) -> bool {
    let name = fixed::record(line);
    for keyword in ["ATOM", "HETATM"] {
        if let Some(rest) = name.strip_prefix(keyword) {
            let rest = rest.trim_start();
            if !rest.is_empty() && rest.bytes().all(|byte| byte.is_ascii_digit()) {
                return true;
            }
        }
    }
    false
}

/// A residue number as an optional 32-bit value.
///
/// The format's sequence field is wider than the sequence columns can hold, so a
/// number outside `i32` is absent rather than wrapped.
#[must_use]
pub(crate) fn sequence_identifier(sequence: i64) -> OptionalI32 {
    match i32::try_from(sequence) {
        Ok(sequence) => OptionalI32::some(sequence),
        Err(_) => OptionalI32::NONE,
    }
}

/// The file-local atom label, or the reserved unlabeled value.
///
/// The serial is a five-column hybrid-36 field. A spelling the scheme cannot
/// decode is reported as unlabeled rather than guessed at.
#[must_use]
pub(crate) fn serial_of(line: &Line<'_>) -> u32 {
    match line.numbers {
        Some(numbers) => numbers.serial,
        None => serial_in(line.text),
    }
}

/// [`serial_of`] for the text of a row.
pub(crate) fn serial_in(text: &str) -> u32 {
    let Some(serial) =
        hybrid36::decode(fixed::raw(text, 7, 11), 5).and_then(|serial| u32::try_from(serial).ok())
    else {
        return 0;
    };
    serial
}

#[cfg(test)]
#[path = "rows_tests.rs"]
mod tests;

/// Columns 79-80 as a formal charge: a digit and a sign in either order.
/// Blank means the file does not say, not that the atom is neutral.
pub(in crate::reader) fn charge_of(line: &Line<'_>) -> (i8, Presence) {
    match line.numbers {
        Some(numbers) => numbers.charge,
        None => charge_in(line.text),
    }
}

/// [`charge_of`] for the text of a row.
pub(in crate::reader) fn charge_in(text: &str) -> (i8, Presence) {
    let field = fixed::text(text, 79, 80);
    let mut magnitude = None;
    let mut sign = None;
    for character in field.chars() {
        match character {
            '+' => sign = Some(1i8),
            '-' => sign = Some(-1i8),
            digit => magnitude = digit.to_digit(10).and_then(|d| i8::try_from(d).ok()),
        }
    }
    match (magnitude, sign) {
        (Some(value), Some(sign)) => (value * sign, Presence::Present),
        _ => (0, Presence::Unknown),
    }
}
