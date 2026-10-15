//! Diagnostics for impossible assembly state and representational limits.

use molframe_core::{Code, Diagnostic};

/// Creates a diagnostic for a missing explicit source model number.
pub(super) fn missing_model_number(model: usize) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_context("model", model.to_string())
        .with_context("required", "explicit source model number")
}

/// Creates the internal-invariant diagnostic used for impossible assembly state.
pub(super) fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}

/// Creates a representational-capacity diagnostic for the named resource.
pub(super) fn capacity(name: &str) -> Diagnostic {
    Diagnostic::new(Code::E1901).with_context("limit", name)
}

/// Wraps one diagnostic for APIs that report a collection of findings.
pub(super) fn single(finding: Diagnostic) -> Vec<Diagnostic> {
    vec![finding]
}
