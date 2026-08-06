//! A finding about data, and the list they accumulate into.
//!
//! A finding always carries a code, and through the code a cause and a remedy —
//! neither is something a call site can forget. The message is optional and
//! defaults to the code's registered cause, so raising a finding that needs no
//! elaboration allocates nothing.

use super::code::{Code, Severity, Strictness};
use crate::span::ByteSpan;
use std::fmt;

/// A named piece of supporting detail attached to a finding.
///
/// The label is fixed by the site that raises the finding; the value comes from
/// the data. Together they carry the specifics a static remedy cannot.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ContextItem {
    label: &'static str,
    value: Box<str>,
}

impl ContextItem {
    /// Creates a context item.
    #[must_use]
    pub fn new(label: &'static str, value: impl Into<Box<str>>) -> Self {
        Self {
            label,
            value: value.into(),
        }
    }

    /// The fixed label.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        self.label
    }

    /// The value taken from the data.
    #[must_use]
    pub fn value(&self) -> &str {
        &self.value
    }
}

/// One finding about the data.
///
/// A finding always has a code, and through the code always has a cause and a
/// remedy — neither is something a call site can forget to supply. The message
/// is optional and defaults to the code's registered cause, so raising a finding
/// that needs no elaboration allocates nothing.
///
/// # Examples
///
/// ```
/// use pdbiox_core::{Code, Diagnostic, Severity};
///
/// let finding = Diagnostic::new(Code::W3203).in_category("atom_site").at_row(17);
/// assert_eq!(finding.severity(), Severity::Loose);
/// assert_eq!(finding.message(), "element inferred from the atom name");
/// assert!(!finding.remedy().is_empty());
/// ```
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Diagnostic {
    code: Code,
    severity: Severity,
    message: Option<Box<str>>,
    span: Option<ByteSpan>,
    category: Option<Box<str>>,
    field: Option<Box<str>>,
    row: Option<u32>,
    context: Vec<ContextItem>,
}

impl Diagnostic {
    /// Creates a finding carrying the code's registered severity and cause.
    #[must_use]
    pub fn new(code: Code) -> Self {
        Self {
            code,
            severity: code.severity(),
            message: None,
            span: None,
            category: None,
            field: None,
            row: None,
            context: Vec::new(),
        }
    }

    /// Replaces the message with one describing this occurrence.
    ///
    /// Use it when the specifics matter; leave it alone when the registered
    /// cause already says everything there is to say.
    #[must_use]
    pub fn with_message(mut self, message: impl Into<Box<str>>) -> Self {
        self.message = Some(message.into());
        self
    }

    /// Overrides the severity for this occurrence.
    ///
    /// The registry decides severity in general. This exists for the cases where
    /// the same code is genuinely worse in one context than another.
    #[must_use]
    pub const fn with_severity(mut self, severity: Severity) -> Self {
        self.severity = severity;
        self
    }

    /// Records where in the source this finding was raised.
    #[must_use]
    pub const fn at(mut self, span: ByteSpan) -> Self {
        self.span = Some(span);
        self
    }

    /// Records the category this finding concerns.
    #[must_use]
    pub fn in_category(mut self, category: impl Into<Box<str>>) -> Self {
        self.category = Some(category.into());
        self
    }

    /// Records the item this finding concerns.
    #[must_use]
    pub fn about_field(mut self, field: impl Into<Box<str>>) -> Self {
        self.field = Some(field.into());
        self
    }

    /// Records the row this finding concerns.
    #[must_use]
    pub const fn at_row(mut self, row: u32) -> Self {
        self.row = Some(row);
        self
    }

    /// Attaches a labelled detail.
    #[must_use]
    pub fn with_context(mut self, label: &'static str, value: impl Into<Box<str>>) -> Self {
        self.context.push(ContextItem::new(label, value));
        self
    }

    /// The stable code.
    #[must_use]
    pub const fn code(&self) -> Code {
        self.code
    }

    /// How bad this finding is.
    #[must_use]
    pub const fn severity(&self) -> Severity {
        self.severity
    }

    /// What happened.
    #[must_use]
    pub fn message(&self) -> &str {
        match &self.message {
            Some(message) => message,
            None => self.code.cause(),
        }
    }

    /// What to do about it.
    #[must_use]
    pub fn remedy(&self) -> &'static str {
        self.code.remedy()
    }

    /// Where in the source this was raised, if a position was known.
    #[must_use]
    pub const fn span(&self) -> Option<ByteSpan> {
        self.span
    }

    /// The category this concerns, if any.
    #[must_use]
    pub fn category(&self) -> Option<&str> {
        self.category.as_deref()
    }

    /// The item this concerns, if any.
    #[must_use]
    pub fn field(&self) -> Option<&str> {
        self.field.as_deref()
    }

    /// The row this concerns, if any.
    #[must_use]
    pub const fn row(&self) -> Option<u32> {
        self.row
    }

    /// The labelled details attached to this finding.
    #[must_use]
    pub fn context(&self) -> &[ContextItem] {
        &self.context
    }

    /// Returns true when this finding is an error at the given strictness.
    #[must_use]
    pub const fn is_error(&self, strictness: Strictness) -> bool {
        self.severity.is_error(strictness)
    }

    /// The key findings are ordered by: worst first, then earliest in the file,
    /// then by code.
    fn order_key(&self) -> (std::cmp::Reverse<Severity>, u32, Code) {
        let offset = self.span.map_or(u32::MAX, |span| span.start.byte_offset);
        (std::cmp::Reverse(self.severity), offset, self.code)
    }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}[{}]: {}",
            self.severity.label(),
            self.code,
            self.message()
        )
    }
}

/// A growing list of findings that keeps itself in reproducible order.
///
/// Callers push in whatever order problems are noticed; [`Diagnostics::finish`]
/// imposes the order that makes two runs comparable.
#[derive(Clone, Default, Debug)]
pub struct Diagnostics {
    findings: Vec<Diagnostic>,
    worst: Option<Severity>,
}

impl Diagnostics {
    /// Creates an empty list.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            findings: Vec::new(),
            worst: None,
        }
    }

    /// Creates a list with room for `capacity` findings.
    ///
    /// Worth doing when a permissive parse over a known-messy file is expected
    /// to raise many, so the list does not grow by repeated reallocation.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            findings: Vec::with_capacity(capacity),
            worst: None,
        }
    }

    /// Records a finding.
    pub fn push(&mut self, finding: Diagnostic) {
        self.worst = Some(match self.worst {
            Some(worst) if worst >= finding.severity => worst,
            _ => finding.severity,
        });
        self.findings.push(finding);
    }

    /// The worst severity recorded so far, if anything has been.
    #[must_use]
    pub const fn worst(&self) -> Option<Severity> {
        self.worst
    }

    /// Returns true when any recorded finding is an error at this strictness.
    #[must_use]
    pub fn has_error(&self, strictness: Strictness) -> bool {
        self.worst.is_some_and(|worst| worst.is_error(strictness))
    }

    /// The findings recorded so far, in the order they were pushed.
    #[must_use]
    pub fn as_slice(&self) -> &[Diagnostic] {
        &self.findings
    }

    /// Returns true when nothing has been recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.findings.is_empty()
    }

    /// The number of findings recorded.
    #[must_use]
    pub fn len(&self) -> usize {
        self.findings.len()
    }

    /// Consumes the list, returning the findings in reproducible order.
    ///
    /// The sort is stable, so two findings identical on every ordering key stay
    /// in the order they were raised.
    #[must_use]
    pub fn finish(mut self) -> Vec<Diagnostic> {
        self.findings.sort_by_key(Diagnostic::order_key);
        self.findings
    }
}

impl Extend<Diagnostic> for Diagnostics {
    fn extend<T: IntoIterator<Item = Diagnostic>>(&mut self, iter: T) {
        for finding in iter {
            self.push(finding);
        }
    }
}

#[cfg(test)]
#[path = "finding_tests.rs"]
mod tests;
