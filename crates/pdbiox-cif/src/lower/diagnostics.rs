//! Source-location helpers for lowering diagnostics.

use pdbiox_core::diagnostic::Diagnostic;

/// Attaches a row when its index fits the public diagnostic representation.
///
/// Oversized input remains diagnosable without truncating its location.
pub(super) fn at_source_row(diagnostic: Diagnostic, row: usize) -> Diagnostic {
    match u32::try_from(row) {
        Ok(row) => diagnostic.at_row(row),
        Err(_) => diagnostic.with_context("source row", row.to_string()),
    }
}
