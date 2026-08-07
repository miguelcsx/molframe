//! Printing results and findings.
//!
//! Results go to standard output and findings to standard error, without
//! exception, so redirecting the one never picks up the other. Findings are
//! rendered by the library rather than reformatted here, because there is one
//! renderer and an error should read the same wherever it surfaces.

use pdbiox::{Diagnostic, ParseMode, Rendered};
use std::fmt::Write as _;
use std::io::Write as _;

/// Machine or human result representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputKind {
    /// Human-readable prose.
    Text,
    /// One structured JSON object.
    Json,
    /// Comma-separated rows.
    Csv,
    /// Tab-separated rows.
    Tsv,
}

/// What the caller asked for that every command needs to know.
#[derive(Clone, Copy, Debug)]
pub struct Context {
    /// Print results as structured data rather than prose.
    pub format: OutputKind,
    /// Suppress findings.
    pub quiet: bool,
    /// Colour the findings.
    pub color: bool,
    /// How much irregularity a read tolerates.
    pub mode: ParseMode,
}

impl Context {
    /// Whether structured JSON was requested.
    #[must_use]
    pub const fn is_json(self) -> bool {
        matches!(self.format, OutputKind::Json)
    }

    /// Delimiter for a requested tabular format.
    #[must_use]
    pub const fn delimiter(self) -> Option<char> {
        match self.format {
            OutputKind::Csv => Some(','),
            OutputKind::Tsv => Some('\t'),
            OutputKind::Text | OutputKind::Json => None,
        }
    }

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

/// A deterministic delimited table.
#[derive(Debug)]
pub struct Table {
    delimiter: char,
    output: String,
}

impl Table {
    /// Starts a table and writes its header.
    #[must_use]
    pub fn new(delimiter: char, header: &[&str]) -> Self {
        let mut table = Self {
            delimiter,
            output: String::new(),
        };
        table.row(header.iter().copied());
        table
    }

    /// Appends one row, quoting fields according to the selected delimiter.
    pub fn row<'a>(&mut self, values: impl IntoIterator<Item = &'a str>) {
        for (position, value) in values.into_iter().enumerate() {
            if position != 0 {
                self.output.push(self.delimiter);
            }
            delimited_field(&mut self.output, value, self.delimiter);
        }
        self.output.push('\n');
    }

    /// Returns the complete table without an extra trailing newline.
    #[must_use]
    pub fn finish(mut self) -> String {
        if self.output.ends_with('\n') {
            self.output.pop();
        }
        self.output
    }
}

fn delimited_field(output: &mut String, value: &str, delimiter: char) {
    if !value.contains([delimiter, '"', '\n', '\r']) {
        output.push_str(value);
        return;
    }
    output.push('"');
    for character in value.chars() {
        if character == '"' {
            output.push_str("\"\"");
        } else {
            output.push(character);
        }
    }
    output.push('"');
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
