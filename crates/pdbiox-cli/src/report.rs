//! Printing results and findings.
//!
//! Results go to standard output and findings to standard error, without
//! exception, so redirecting the one never picks up the other. Findings are
//! rendered by the library rather than reformatted here, because there is one
//! renderer and an error should read the same wherever it surfaces.

use pdbiox::{Diagnostic, ParseMode, Rendered};
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

mod rows;

pub use rows::RowWriter;

/// Machine or human result representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputKind {
    /// Human-readable prose.
    Text,
    /// One structured JSON object.
    Json,
    /// One structured object per line for streaming.
    JsonLines,
    /// Comma-separated rows.
    Csv,
    /// Tab-separated rows.
    Tsv,
    /// Apache Arrow IPC file output.
    Arrow,
    /// Apache Parquet file output.
    Parquet,
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
    /// Effective analysis policy after config and CLI overrides.
    pub policy: &'static pdbiox::AnalysisPolicy,
    /// Optional destination for the primary rendered result.
    pub output: Option<&'static Path>,
    /// Optional destination for a standalone provenance record.
    pub provenance: Option<&'static Path>,
    /// Configured Chemical Component Dictionary, if any.
    pub ccd: Option<&'static Path>,
    /// Exact release identifier paired with the configured dictionary.
    pub ccd_version: Option<&'static str>,
    /// Shared native executor, memory, cancellation, scratch and spill policy.
    pub execution: &'static pdbiox::core::ExecutionContext,
    /// Explicit missing-element behavior used by every structural reader.
    pub missing_element_policy: pdbiox::MissingElementPolicy,
    /// Explicit behavior for identifiers that cannot delimit adjacent residues.
    pub residue_boundary_policy: pdbiox::AmbiguousResidueBoundaryPolicy,
}

impl Context {
    /// Verifies result destinations before scientific work begins.
    pub fn prepare_outputs(self) -> std::io::Result<()> {
        for path in [self.output, self.provenance].into_iter().flatten() {
            let directory = path
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty());
            let _temporary = match directory {
                Some(parent) => tempfile::NamedTempFile::new_in(parent)?,
                None => tempfile::NamedTempFile::new_in(".")?,
            };
        }
        Ok(())
    }

    /// Whether structured JSON was requested.
    #[must_use]
    pub const fn is_json(self) -> bool {
        matches!(self.format, OutputKind::Json | OutputKind::JsonLines)
    }

    /// Delimiter for a requested tabular format.
    #[must_use]
    pub const fn delimiter(self) -> Option<char> {
        match self.format {
            OutputKind::Csv => Some(','),
            OutputKind::Tsv => Some('\t'),
            OutputKind::Text
            | OutputKind::Json
            | OutputKind::JsonLines
            | OutputKind::Arrow
            | OutputKind::Parquet => None,
        }
    }

    /// Whether a binary Arrow-backed file was requested.
    #[must_use]
    pub const fn is_table_file(self) -> bool {
        matches!(self.format, OutputKind::Arrow | OutputKind::Parquet)
    }

    /// Renders a sequence of JSON objects as one array or a JSON Lines stream.
    #[must_use]
    pub fn json_records(self, records: &[String]) -> String {
        match self.format {
            OutputKind::JsonLines => records.join("\n"),
            _ => json_array(records),
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
    pub fn result(self, text: &str) {
        let rendered = self.with_embedded_provenance(text);
        if let Some(path) = self.output {
            if let Err(error) = std::fs::write(path, format!("{rendered}\n")) {
                eprintln!("could not write result to {}: {error}", path.display());
            }
        } else {
            let mut out = std::io::stdout().lock();
            let _ = writeln!(out, "{rendered}");
        }
        if let Some(path) = self.provenance
            && let Err(error) = std::fs::write(path, format!("{}\n", self.provenance_json()))
        {
            eprintln!("could not write provenance to {}: {error}", path.display());
        }
    }

    fn with_embedded_provenance(self, text: &str) -> String {
        let provenance = self.provenance_json();
        match self.format {
            OutputKind::Json if text.starts_with('{') && text.ends_with('}') => {
                let Some(body) = text.strip_suffix('}') else {
                    return text.to_owned();
                };
                let separator = if body.len() == 1 { "" } else { "," };
                format!("{body}{separator}\"_provenance\":{provenance}}}")
            }
            OutputKind::Json => format!("{{\"result\":{text},\"_provenance\":{provenance}}}"),
            OutputKind::JsonLines => {
                format!("{text}\n{{\"_provenance\":{provenance}}}")
            }
            OutputKind::Csv | OutputKind::Tsv => {
                format!("{text}\n# pdbiox-provenance={provenance}")
            }
            _ => text.to_owned(),
        }
    }

    pub(super) fn provenance_json(self) -> String {
        let mut record = Json::new();
        record
            .text("pdbiox_version", env!("CARGO_PKG_VERSION"))
            .text(
                "profile",
                &self
                    .policy
                    .profile()
                    .map_or_else(|| "modified".to_owned(), |profile| profile.to_string()),
            )
            .text("policy_fingerprint", &self.policy.fingerprint().to_string())
            .text("policy", &self.policy.to_string())
            .text(
                "invocation",
                &std::env::args().collect::<Vec<_>>().join(" "),
            )
            .number("workers", self.execution.worker_budget())
            .number("memory_budget", self.execution.memory_budget().bytes())
            .number(
                "spill_budget",
                self.execution.temp_storage_policy().max_bytes(),
            );
        record.text(
            "missing_element_policy",
            match self.missing_element_policy {
                pdbiox::MissingElementPolicy::PreserveUnknown => "preserve-unknown",
                pdbiox::MissingElementPolicy::InferFromAtomName => "infer-from-atom-name",
            },
        );
        record.text(
            "ambiguous_residue_boundary_policy",
            match self.residue_boundary_policy {
                pdbiox::AmbiguousResidueBoundaryPolicy::Reject => "reject",
                pdbiox::AmbiguousResidueBoundaryPolicy::InferFromFileOrder => {
                    "infer-from-file-order"
                }
            },
        );
        if let Some(path) = self.ccd {
            record.text("ccd", &path.display().to_string());
        }
        if let Some(version) = self.ccd_version {
            record.text("ccd_version", version);
        }
        record.finish()
    }

    /// Provenance fields suitable for Arrow schema metadata.
    #[must_use]
    pub fn provenance_metadata(self) -> std::collections::BTreeMap<String, String> {
        let mut metadata = std::collections::BTreeMap::from([
            (
                "pdbiox.version".to_owned(),
                env!("CARGO_PKG_VERSION").to_owned(),
            ),
            (
                "pdbiox.policy_fingerprint".to_owned(),
                self.policy.fingerprint().to_string(),
            ),
            ("pdbiox.policy".to_owned(), self.policy.to_string()),
            (
                "pdbiox.missing_element_policy".to_owned(),
                match self.missing_element_policy {
                    pdbiox::MissingElementPolicy::PreserveUnknown => "preserve-unknown",
                    pdbiox::MissingElementPolicy::InferFromAtomName => "infer-from-atom-name",
                }
                .to_owned(),
            ),
            (
                "pdbiox.ambiguous_residue_boundary_policy".to_owned(),
                match self.residue_boundary_policy {
                    pdbiox::AmbiguousResidueBoundaryPolicy::Reject => "reject",
                    pdbiox::AmbiguousResidueBoundaryPolicy::InferFromFileOrder => {
                        "infer-from-file-order"
                    }
                }
                .to_owned(),
            ),
        ]);
        if let Some(path) = self.ccd {
            metadata.insert("pdbiox.ccd".to_owned(), path.display().to_string());
        }
        if let Some(version) = self.ccd_version {
            metadata.insert("pdbiox.ccd_version".to_owned(), version.to_owned());
        }
        metadata
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

pub(super) fn delimited_field(output: &mut String, value: &str, delimiter: char) {
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
pub(super) fn escape(text: &str) -> String {
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
