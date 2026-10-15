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
//! A binary table carries the same header: its magic, a little-endian `u16`
//! header length, the header text, then the payload.
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

/// A sink that accepts the identity headers and nothing else.
#[derive(Clone, Copy, Debug, Default)]
pub struct IdentityOnly;

impl TableLines for IdentityOnly {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, _: &str) -> Result<(), TableTextError> {
        Err(TableTextError::at(
            "table header",
            line,
            "only headers are allowed here",
        ))
    }
}

/// Splits a binary table into its declared identity and its payload, handing any
/// other header of the kind to `lines`.
///
/// # Errors
///
/// [`TableTextError`] for a wrong magic, a truncated header or a bad identity.
pub fn split_binary<'a, L: TableLines>(
    what: &'static str,
    magic: [u8; 4],
    bytes: &'a [u8],
    lines: &mut L,
) -> Result<(TableIdentity, &'a [u8]), TableTextError> {
    let whole = |message: &str| TableTextError::at(what, 0, message);
    if bytes.get(..4) != Some(&magic) {
        return Err(whole("not a table of this kind"));
    }
    let (Some(&low), Some(&high)) = (bytes.get(4), bytes.get(5)) else {
        return Err(whole("truncated before the header length"));
    };
    let length = usize::from(u16::from_le_bytes([low, high]));
    let Some(header) = bytes.get(6..6 + length) else {
        return Err(whole("truncated inside the header"));
    };
    let identity = read_table(what, text_of(what, header)?, lines)?;
    let payload = match bytes.get(6 + length..) {
        Some(payload) => payload,
        None => &[],
    };
    Ok((identity, payload))
}

/// The inverse of [`split_binary`]: `payload` behind a header declaring `identity`
/// followed by `extra` header lines.
///
/// Returns `None` when the header text does not fit the `u16` length.
#[must_use]
pub fn join_binary(
    magic: [u8; 4],
    identity: &TableIdentity,
    extra: &str,
    payload: &[u8],
) -> Option<Vec<u8>> {
    let header = format!(
        "name: {}\nversion: {}\nsource: {}\nlicence: {}\n{extra}",
        identity.name, identity.version, identity.source, identity.licence
    );
    let length = u16::try_from(header.len()).ok()?;
    let mut out = Vec::with_capacity(6 + header.len() + payload.len());
    out.extend_from_slice(&magic);
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    out.extend_from_slice(payload);
    Some(out)
}

#[cfg(test)]
#[path = "table_text_tests.rs"]
mod tests;
