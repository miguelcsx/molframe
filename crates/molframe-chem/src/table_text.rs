//! The one text format every reference table shares, and its reader.
//!
//! A table file is a versioned header followed by rows:
//!
//! ```text
//! # SPDX-License-Identifier: CC0-1.0
//! format: 1
//! name: lab-radii
//! version: 2
//! source: Smith lab, internal survey
//! licence: CC0-1.0
//!
//! # rows, in the layout of the table's own kind
//! C 1.75
//! ```
//!
//! `name`, `version`, `source` and `licence` become the [`TableIdentity`]; `format`
//! is optional and, when present, must be `1`. Any other `key: value` line is
//! offered to the table kind (a base set, a residue class), everything else is a
//! row. Blank lines and `#` comments are skipped. A problem is reported with its
//! one-based line number and a table is accepted whole or not at all.
//!
//! Cost: one pass over the text, no allocation beyond the four identity strings
//! and whatever the kind keeps from its rows.

use crate::TableIdentity;
use core::fmt;

/// The version of the text format this reader understands.
pub const TABLE_FORMAT_VERSION: u32 = 1;

/// Why a table was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TableTextError {
    what: &'static str,
    line: usize,
    message: String,
}

impl TableTextError {
    /// A problem on one-based `line` of a `what` table; zero for the whole table.
    pub fn at(what: &'static str, line: usize, message: impl Into<String>) -> Self {
        Self {
            what,
            line,
            message: message.into(),
        }
    }

    /// The one-based line the problem is on, or zero for a whole-table problem.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

impl fmt::Display for TableTextError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(out, "{}: {}", self.what, self.message)
        } else {
            write!(out, "{}, line {}: {}", self.what, self.line, self.message)
        }
    }
}

impl core::error::Error for TableTextError {}

/// What a table kind does with the lines the shared reader does not own.
pub trait TableLines {
    /// Reads a `key: value` header of this kind; `Ok(false)` when it is not one.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when the header is one of this kind's but its value is
    /// refused.
    fn header(&mut self, line: usize, key: &str, value: &str) -> Result<bool, TableTextError>;

    /// Reads one row.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when the row is malformed.
    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError>;
}

/// Reads `text`, handing headers and rows of the kind to `lines`, and returns the
/// declared identity.
///
/// # Errors
///
/// [`TableTextError`] for a missing, repeated or unknown header, an unsupported
/// `format`, or whatever `lines` refuses.
pub fn read_table<L: TableLines>(
    what: &'static str,
    text: &str,
    lines: &mut L,
) -> Result<TableIdentity, TableTextError> {
    let mut fields: [Option<&str>; 4] = [None; 4];
    let mut format_seen = false;
    for (index, raw) in text.lines().enumerate() {
        let number = index + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((key, value)) = line.split_once(':') else {
            lines.row(number, line)?;
            continue;
        };
        let (key, value) = (key.trim(), value.trim());
        if key == "format" {
            if format_seen || value.parse::<u32>() != Ok(TABLE_FORMAT_VERSION) {
                return Err(TableTextError::at(
                    what,
                    number,
                    format!("`format` must be declared once, as {TABLE_FORMAT_VERSION}"),
                ));
            }
            format_seen = true;
            continue;
        }
        let slot = match key {
            "name" => Some(0),
            "version" => Some(1),
            "source" => Some(2),
            "licence" | "license" => Some(3),
            _ => None,
        };
        let Some(slot) = slot else {
            if !lines.header(number, key, value)? {
                return Err(TableTextError::at(
                    what,
                    number,
                    format!("unknown header `{key}`"),
                ));
            }
            continue;
        };
        if value.is_empty() {
            return Err(TableTextError::at(
                what,
                number,
                format!("`{key}` has no value"),
            ));
        }
        match fields.get_mut(slot) {
            Some(field @ None) => *field = Some(value),
            _ => {
                return Err(TableTextError::at(
                    what,
                    number,
                    format!("`{key}` is declared twice"),
                ));
            }
        }
    }
    let [name, version, source, licence] = fields;
    Ok(TableIdentity::supplied(
        need(what, name, "name")?,
        need(what, version, "version")?,
        need(what, source, "source")?,
        need(what, licence, "licence")?,
    ))
}

fn need<'a>(
    what: &'static str,
    value: Option<&'a str>,
    key: &str,
) -> Result<&'a str, TableTextError> {
    value.ok_or_else(|| TableTextError::at(what, 0, format!("missing `{key}` header")))
}

/// The text of a table file, or why it is not text.
///
/// # Errors
///
/// [`TableTextError`] when `bytes` are not UTF-8.
pub fn text_of<'a>(what: &'static str, bytes: &'a [u8]) -> Result<&'a str, TableTextError> {
    core::str::from_utf8(bytes).map_err(|error| TableTextError::at(what, 0, error.to_string()))
}

#[cfg(test)]
#[path = "table_text_tests.rs"]
mod tests;
