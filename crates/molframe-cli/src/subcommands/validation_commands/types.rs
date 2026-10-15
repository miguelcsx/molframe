//! The shapes shared by every validation check: its request, subject and rows.

use crate::args::ValidateArguments;
use crate::report::Context;

/// What `validate` was asked to check, and with which explicit parameters.
pub(crate) type ValidationOptions<'a> = &'a ValidateArguments;

/// The inputs every check can draw on.
pub(super) struct Subject<'a> {
    pub(super) structure: &'a molframe::Structure,
    pub(super) options: ValidationOptions<'a>,
    pub(super) context: Context,
}

#[derive(Debug)]
pub(super) struct Row {
    pub(super) check: String,
    pub(super) item: String,
    pub(super) value: String,
    /// Whether the row reports a problem; summaries and provenance do not.
    pub(super) finding: bool,
}

impl Row {
    pub(super) fn new(
        check: impl Into<String>,
        item: impl Into<String>,
        value: impl Into<String>,
    ) -> Self {
        Self {
            check: check.into(),
            item: item.into(),
            value: value.into(),
            finding: true,
        }
    }

    /// A row that informs without counting as a finding.
    pub(super) fn info(
        check: impl Into<String>,
        item: impl Into<String>,
        value: impl Into<String>,
    ) -> Self {
        Self {
            finding: false,
            ..Self::new(check, item, value)
        }
    }
}
