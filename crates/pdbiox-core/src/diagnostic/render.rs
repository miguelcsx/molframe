//! Presentation of findings.
//!
//! There is exactly one renderer. The command line, a Rust `Display` and any
//! binding all route through it, so a finding reads the same wherever it
//! surfaces and a user searching for the text of an error finds one answer.
//!
//! Rendering borrows everything it needs and writes straight to the formatter,
//! so showing a finding costs no allocation even when it quotes the source line.

use super::{Diagnostic, Severity};
use std::fmt;

const RESET: &str = "\u{1b}[0m";
const BOLD: &str = "\u{1b}[1m";
const RED: &str = "\u{1b}[31m";
const YELLOW: &str = "\u{1b}[33m";
const BLUE: &str = "\u{1b}[34m";

/// A finding prepared for display.
///
/// Built with [`Diagnostic::render`](Diagnostic) style helpers on this type and
/// then formatted. Without a source buffer it renders the heading and the
/// labelled details; with one it also quotes and underlines the offending line.
///
/// # Examples
///
/// ```
/// use pdbiox_core::{Code, Diagnostic, Rendered};
///
/// let finding = Diagnostic::new(Code::E1106);
/// let text = Rendered::new(&finding).to_string();
/// assert!(text.starts_with("error[PDBIOX-E1106]"));
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Rendered<'a> {
    finding: &'a Diagnostic,
    source: Option<&'a [u8]>,
    origin: Option<&'a str>,
    color: bool,
}

impl<'a> Rendered<'a> {
    /// Prepares a finding for display with no source context and no colour.
    #[must_use]
    pub const fn new(finding: &'a Diagnostic) -> Self {
        Self {
            finding,
            source: None,
            origin: None,
            color: false,
        }
    }

    /// Supplies the bytes the finding was raised against, so the offending line
    /// can be quoted.
    #[must_use]
    pub const fn with_source(mut self, source: &'a [u8]) -> Self {
        self.source = Some(source);
        self
    }

    /// Names where the input came from, for the location line.
    #[must_use]
    pub const fn with_origin(mut self, origin: &'a str) -> Self {
        self.origin = Some(origin);
        self
    }

    /// Enables ANSI colour.
    ///
    /// Callers decide this from whether the stream is a terminal and from the
    /// user's own preference; the renderer does not probe the environment.
    #[must_use]
    pub const fn with_color(mut self, color: bool) -> Self {
        self.color = color;
        self
    }

    fn accent(self) -> &'static str {
        if !self.color {
            return "";
        }
        match self.finding.severity() {
            Severity::Invalidating | Severity::Breaking => RED,
            Severity::Strict | Severity::Loose => YELLOW,
            Severity::Info => BLUE,
        }
    }

    fn paint(self, code: &'static str) -> &'static str {
        if self.color { code } else { "" }
    }

    /// Writes the `error[CODE]: message` heading.
    fn write_heading(self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let finding = self.finding;
        write!(
            f,
            "{}{}{}[{}]{}: {}{}{}",
            self.accent(),
            self.paint(BOLD),
            finding.severity().label(),
            finding.code(),
            self.paint(RESET),
            self.paint(BOLD),
            finding.message(),
            self.paint(RESET),
        )
    }

    /// Writes the ` ┌─ origin:line:column` locator.
    fn write_location(self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (Some(span), gutter) = (self.finding.span(), self.paint(BLUE)) else {
            if let Some(origin) = self.origin {
                write!(
                    f,
                    "\n  {}┌─{} {origin}",
                    self.paint(BLUE),
                    self.paint(RESET)
                )?;
            }
            return Ok(());
        };
        write!(f, "\n  {gutter}┌─{} ", self.paint(RESET))?;
        match self.origin {
            Some(origin) => write!(f, "{origin}:{}", span.start),
            None => write!(f, "{}", span.start),
        }
    }

    /// Quotes the source line the span falls on and underlines the span.
    fn write_snippet(self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (Some(source), Some(span)) = (self.source, self.finding.span()) else {
            return Ok(());
        };
        let Some(line) = line_containing(source, span.start.byte_offset as usize) else {
            return Ok(());
        };
        let Ok(text) = str::from_utf8(line.text) else {
            return Ok(());
        };

        let gutter = self.paint(BLUE);
        let reset = self.paint(RESET);
        let number = span.start.line;
        let width = decimal_width(number);

        write!(f, "\n {:width$} {gutter}│{reset}", "")?;
        write!(f, "\n {number} {gutter}│{reset} {}", text.trim_end())?;

        let column = (span.start.byte_offset as usize).saturating_sub(line.start);
        let caret_count = (span.len() as usize).clamp(1, text.len().saturating_sub(column).max(1));
        write!(f, "\n {:width$} {gutter}│{reset} {:column$}", "", "")?;
        write!(f, "{}", self.accent())?;
        for _ in 0..caret_count {
            f.write_str("^")?;
        }
        f.write_str(reset)
    }

    /// Writes the labelled details and the remedy.
    fn write_details(self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let finding = self.finding;
        let gutter = self.paint(BLUE);
        let reset = self.paint(RESET);

        if let Some(category) = finding.category() {
            write!(f, "\n  {gutter}={reset} category: {category}")?;
        }
        if let Some(field) = finding.field() {
            write!(f, "\n  {gutter}={reset} item: {field}")?;
        }
        if let Some(row) = finding.row() {
            write!(f, "\n  {gutter}={reset} row: {row}")?;
        }
        for item in finding.context() {
            write!(f, "\n  {gutter}={reset} {}: {}", item.label(), item.value())?;
        }
        write!(f, "\n  {gutter}={reset} help: {}", finding.remedy())
    }
}

impl fmt::Display for Rendered<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write_heading(f)?;
        self.write_location(f)?;
        self.write_snippet(f)?;
        self.write_details(f)
    }
}

/// The bounds of the line containing a byte offset.
struct Line<'a> {
    text: &'a [u8],
    start: usize,
}

fn line_containing(source: &[u8], offset: usize) -> Option<Line<'_>> {
    if offset > source.len() {
        return None;
    }
    let start = source[..offset]
        .iter()
        .rposition(|&b| b == b'\n')
        .map_or(0, |at| at + 1);
    let end = source[start..]
        .iter()
        .position(|&b| b == b'\n')
        .map_or(source.len(), |at| start + at);
    Some(Line {
        text: source.get(start..end)?,
        start,
    })
}

const fn decimal_width(mut value: u32) -> usize {
    let mut width = 1;
    while value >= 10 {
        value /= 10;
        width += 1;
    }
    width
}

#[cfg(test)]
#[path = "render_tests.rs"]
mod tests;
