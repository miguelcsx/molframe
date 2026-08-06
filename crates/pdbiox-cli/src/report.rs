//! Printing results and findings.
//!
//! Results go to standard output and findings to standard error, without
//! exception, so redirecting the one never picks up the other. Findings are
//! rendered by the library rather than reformatted here, because there is one
//! renderer and an error should read the same wherever it surfaces.

use pdbiox::{Diagnostic, ParseMode, Rendered};
use std::fmt::Write as _;
use std::io::Write as _;

/// What the caller asked for that every command needs to know.
#[derive(Clone, Copy, Debug)]
pub struct Context {
    /// Print results as structured data rather than prose.
    pub format: bool,
    /// Suppress findings.
    pub quiet: bool,
    /// Colour the findings.
    pub color: bool,
    /// How much irregularity a read tolerates.
    pub mode: ParseMode,
}

impl Context {
    /// Prints findings to standard error.
    pub fn findings(self, findings: &[Diagnostic], origin: &str) {
        if self.quiet {
            return;
        }
        let mut err = std::io::stderr().lock();
        for finding in findings {
            let rendered = Rendered::new(finding)
                .with_origin(origin)
                .with_color(self.color);
            let _ = writeln!(err, "{rendered}");
        }
    }

    /// Prints a result to standard output.
    ///
    /// Takes the context by value so that every printing path goes through it,
    /// even the ones that do not currently consult a field.
    #[allow(clippy::unused_self, reason = "one printing path for every command")]
    pub fn result(self, text: &str) {
        let mut out = std::io::stdout().lock();
        let _ = writeln!(out, "{text}");
    }
}

/// Builds a structured object, one field at a time.
///
/// A dependency on a serialisation library would be the fifth crate in a layer
/// that is meant to have none; the objects this prints have a handful of fields
/// and no nesting to speak of.
#[derive(Debug, Default)]
pub struct Json {
    body: String,
}

impl Json {
    /// Starts an empty object.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a string field.
    pub fn text(&mut self, key: &str, value: &str) -> &mut Self {
        self.separate();
        let _ = write!(self.body, "\"{}\":\"{}\"", escape(key), escape(value));
        self
    }

    /// Adds a numeric field.
    pub fn number(&mut self, key: &str, value: impl std::fmt::Display) -> &mut Self {
        self.separate();
        let _ = write!(self.body, "\"{}\":{value}", escape(key));
        self
    }

    /// Adds a field holding an already-built list.
    pub fn raw(&mut self, key: &str, value: &str) -> &mut Self {
        self.separate();
        let _ = write!(self.body, "\"{}\":{value}", escape(key));
        self
    }

    fn separate(&mut self) {
        if !self.body.is_empty() {
            self.body.push(',');
        }
    }

    /// Renders the object.
    #[must_use]
    pub fn finish(&self) -> String {
        format!("{{{}}}", self.body)
    }
}

/// Renders a list of already-built objects.
#[must_use]
pub fn json_array(items: &[String]) -> String {
    format!("[{}]", items.join(","))
}

/// Escapes the characters a string field cannot carry literally.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if control < ' ' => {
                let _ = write!(out, "\\u{:04x}", control as u32);
            }
            other => out.push(other),
        }
    }
    out
}

#[cfg(test)]
#[path = "report_tests.rs"]
mod tests;
