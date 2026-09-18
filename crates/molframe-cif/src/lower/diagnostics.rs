//! Source-location helpers for lowering diagnostics.

use molframe_core::diagnostic::Diagnostic;

/// Attaches the exact logical row without narrowing it to a chunk-local index.
pub(super) fn at_source_row(diagnostic: Diagnostic, row: usize) -> Diagnostic {
    diagnostic.at_row(row as u64)
}
