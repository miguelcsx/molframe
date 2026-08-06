//! Fixed-width fields.
//!
//! Columns are 1-based and inclusive, as the format's own documentation quotes
//! them, and a line shorter than the field is read as though it were padded with
//! spaces — which is what the Fortran that defined the format did, and what
//! files written by that Fortran assume.
//!
//! Nothing here allocates. A field is a borrowed slice of the line.

/// The text of columns `from` to `to`, trimmed, or the empty string when the
/// line stops short.
#[must_use]
pub fn text(line: &str, from: usize, to: usize) -> &str {
    raw(line, from, to).trim()
}

/// The text of columns `from` to `to` without trimming.
///
/// The leading space of an atom name is the only thing distinguishing a
/// one-letter element from a two-letter one, so some fields must not be trimmed
/// before they are interpreted.
#[must_use]
pub fn raw(line: &str, from: usize, to: usize) -> &str {
    let start = from.saturating_sub(1);
    let end = to.min(line.len());
    if start >= end {
        return "";
    }
    match line.get(start..end) {
        Some(field) => field,
        None => "",
    }
}

/// The record name a line carries.
#[must_use]
pub fn record(line: &str) -> &str {
    text(line, 1, 6)
}

/// An integer read from columns `from` to `to`.
#[must_use]
pub fn integer(line: &str, from: usize, to: usize) -> Option<i64> {
    let field = text(line, from, to);
    if field.is_empty() {
        None
    } else {
        field.parse().ok()
    }
}

/// A real number read from columns `from` to `to`.
#[must_use]
pub fn real(line: &str, from: usize, to: usize) -> Option<f64> {
    let field = text(line, from, to);
    if field.is_empty() {
        None
    } else {
        field.parse().ok()
    }
}

#[cfg(test)]
#[path = "fixed_tests.rs"]
mod tests;
