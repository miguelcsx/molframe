//! Writing and reading back the header of a stored table.
//!
//! The stored form is the shared text format: an SPDX line, the identity
//! headers, a comment recording the pinned digest of the source bytes, then the
//! rows. Reading back is the same reader every table kind uses, with a sink that
//! accepts any row, so a stored file is checked for identity without knowing its
//! kind. Cost: linear in the table's size.

use super::manifest::Entry;
use molframe_chem::{TableLines, TableTextError, read_table};

const DIGEST_COMMENT: &str = "# source-sha256: ";
const WHAT: &str = "stored table";

/// The text of `entry` holding `rows`.
pub(super) fn table_text(entry: &Entry, rows: &[String]) -> String {
    let source = match &entry.source {
        Some(source) => source.as_str(),
        None => entry.url.as_str(),
    };
    let mut text = format!(
        "# SPDX-License-Identifier: {licence}\nformat: 1\nname: {name}\nversion: {version}\n\
         source: {source}\nlicence: {licence}\n{DIGEST_COMMENT}{digest}\n\n",
        licence = entry.licence,
        name = entry.name,
        version = entry.version,
        digest = entry.sha256.to_ascii_lowercase(),
    );
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

/// The rows of a local file to append: every line that is neither blank, a
/// comment nor a `key: value` header.
pub(super) fn bare_rows(text: &str) -> Vec<String> {
    text.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#') && !line.contains(':'))
        .map(str::to_owned)
        .collect()
}

/// A sink that accepts any header of a kind and any row, counting the rows.
#[derive(Default)]
pub(super) struct AnyRows {
    pub(super) rows: usize,
}

impl TableLines for AnyRows {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(true)
    }

    fn row(&mut self, _: usize, _: &str) -> Result<(), TableTextError> {
        self.rows += 1;
        Ok(())
    }
}

/// What a stored table declares.
pub(super) struct Stored {
    pub(super) name: String,
    pub(super) version: String,
    pub(super) licence: String,
    pub(super) digest: Option<String>,
    pub(super) rows: usize,
}

/// Parses `text` with the shared reader.
pub(super) fn read_back(text: &str) -> Result<Stored, TableTextError> {
    let mut sink = AnyRows::default();
    let identity = read_table(WHAT, text, &mut sink)?;
    let digest = text
        .lines()
        .find_map(|line| line.strip_prefix(DIGEST_COMMENT))
        .map(|digest| digest.trim().to_owned());
    Ok(Stored {
        name: identity.name.to_string(),
        version: identity.version.to_string(),
        licence: identity.licence.to_string(),
        digest,
        rows: sink.rows,
    })
}

#[cfg(test)]
#[path = "render_tests.rs"]
mod tests;
