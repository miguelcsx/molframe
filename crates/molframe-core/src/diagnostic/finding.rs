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
/// use molframe_core::{Code, Diagnostic, Severity};
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
    row_present: bool,
    row: u64,
    message: Option<Box<str>>,
    span: Option<Box<ByteSpan>>,
    category: Option<Box<str>>,
    field: Option<Box<str>>,
    context: Vec<ContextItem>,
}

impl Diagnostic {
    /// Heap bytes retained by this occurrence's owned details.
    #[must_use]
    pub fn retained_bytes(&self) -> usize {
        let text = self.message.as_ref().map_or(0, |value| value.len())
            + self.category.as_ref().map_or(0, |value| value.len())
            + self.field.as_ref().map_or(0, |value| value.len())
            + self
                .context
                .iter()
                .map(|item| item.value.len())
                .sum::<usize>();
        text.saturating_add(
            self.context
                .capacity()
                .saturating_mul(std::mem::size_of::<ContextItem>()),
        )
        .saturating_add(
            self.span
                .as_ref()
                .map_or(0, |_| std::mem::size_of::<ByteSpan>()),
        )
    }

    /// Creates a finding carrying the code's registered severity and cause.
    #[must_use]
    pub fn new(code: Code) -> Self {
        Self {
            code,
            severity: code.severity(),
            row_present: false,
            row: 0,
            message: None,
            span: None,
            category: None,
            field: None,
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
    pub fn at(mut self, span: ByteSpan) -> Self {
        self.span = Some(Box::new(span));
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
    pub const fn at_row(mut self, row: u64) -> Self {
        self.row = row;
        self.row_present = true;
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
    pub fn span(&self) -> Option<ByteSpan> {
        self.span.as_deref().copied()
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
    pub const fn row(&self) -> Option<u64> {
        if self.row_present {
            Some(self.row)
        } else {
            None
        }
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
    fn order_key(&self) -> (std::cmp::Reverse<Severity>, u64, Code) {
        let offset = self
            .span
            .as_deref()
            .map_or(u64::MAX, |span| span.start.byte_offset);
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
///
/// The list is bounded. Past its ceiling, findings are counted rather than
/// retained, and [`Diagnostics::finish`] appends one `MOLFRAME-W1901` recording
/// how many were suppressed. The worst severity still reflects every finding
/// pushed, retained or not, so a decision taken on severity is unaffected by
/// the ceiling.
#[derive(Clone, Debug)]
pub struct Diagnostics {
    findings: Vec<Diagnostic>,
    worst: Option<Severity>,
    ceiling: usize,
    suppressed: u64,
}

impl Default for Diagnostics {
    fn default() -> Self {
        Self::new()
    }
}

impl Diagnostics {
    /// How many findings a list retains before it starts counting instead.
    ///
    /// A systematically malformed column raises one finding per row, and a
    /// finding carries a message, a category, a field and a context vector —
    /// roughly ninety bytes and up to four allocations each. At a hundred
    /// million rows that is nine gigabytes of diagnostics describing a single
    /// cause, which exhausts memory reporting a problem rather than reporting
    /// it.
    ///
    /// Ten thousand is far past the point where a reader learns anything new
    /// and far below the point where the list itself is the problem.
    pub const DEFAULT_CEILING: usize = 10_000;

    /// Creates an empty list with the default ceiling.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            findings: Vec::new(),
            worst: None,
            ceiling: Self::DEFAULT_CEILING,
            suppressed: 0,
        }
    }

    /// Creates a list with room for `capacity` findings.
    ///
    /// Worth doing when a permissive parse over a known-messy file is expected
    /// to raise many, so the list does not grow by repeated reallocation. The
    /// ceiling is raised to `capacity` where that is higher than the default,
    /// because a caller who has reserved the room has said what they want.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            findings: Vec::with_capacity(capacity),
            worst: None,
            ceiling: capacity.max(Self::DEFAULT_CEILING),
            suppressed: 0,
        }
    }

    /// Sets how many findings this list retains before counting instead.
    #[must_use]
    pub const fn with_ceiling(mut self, ceiling: usize) -> Self {
        self.ceiling = ceiling;
        self
    }

    /// How many findings this list retains before counting instead.
    #[must_use]
    pub const fn ceiling(&self) -> usize {
        self.ceiling
    }

    /// How many findings were counted rather than retained.
    #[must_use]
    pub const fn suppressed(&self) -> u64 {
        self.suppressed
    }

    /// Records a finding.
    ///
    /// Past the ceiling the finding is counted and dropped. The severity is
    /// recorded either way, so nothing that depends on the worst severity
    /// changes with the ceiling.
    pub fn push(&mut self, finding: Diagnostic) {
        self.worst = Some(match self.worst {
            Some(worst) if worst >= finding.severity => worst,
            _ => finding.severity,
        });
        if self.findings.len() >= self.ceiling {
            self.suppressed = self.suppressed.saturating_add(1);
            return;
        }
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
        if self.suppressed != 0 {
            self.findings.push(
                Diagnostic::new(Code::W1901)
                    .with_context("suppressed", self.suppressed.to_string())
                    .with_context("retained", self.findings.len().to_string()),
            );
        }
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
