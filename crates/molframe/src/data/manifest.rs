//! A manifest: the TOML list of tables a user is willing to fetch.
//!
//! ```toml
//! [[entry]]
//! name = "lab-radii"
//! version = "2"
//! url = "https://example.org/radii.tsv"
//! licence = "CC0-1.0"
//! sha256 = "<64 hex digits>"
//! extension = "radii"
//! converter = "tsv"            # optional: tsv, stereo_chemical_props
//! append = ["extra.rows"]      # optional: local row files, relative to the manifest
//! ```
//!
//! Nothing about any dataset is known to the library; every name, address and
//! digest comes from the file. Cost: one parse of a file of a few entries, then
//! linear validation of each field.

use super::error::DataError;
use crate::document::read_document;
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// One fetchable table.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Entry {
    /// Cache entry name; the file is `<name>.<extension>`.
    pub name: String,
    /// Version the table declares.
    pub version: String,
    /// Where the source bytes live; must be `https`.
    pub url: String,
    /// SPDX expression covering the stored table, appended rows included.
    pub licence: String,
    /// Pinned SHA-256 of the downloaded bytes.
    pub sha256: String,
    /// File extension of the cache entry, naming its kind.
    pub extension: String,
    /// Where the data came from; the URL when omitted.
    #[serde(default)]
    pub source: Option<String>,
    /// Converter id; absent means the download is already a complete table.
    #[serde(default)]
    pub converter: Option<String>,
    /// Local row files concatenated onto the converted table.
    #[serde(default)]
    pub append: Vec<PathBuf>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct File {
    entry: Vec<Entry>,
}

/// A validated manifest and the directory its relative paths resolve against.
#[derive(Clone, Debug)]
pub struct Manifest {
    path: PathBuf,
    entries: Vec<Entry>,
}

impl Manifest {
    /// Reads and validates the manifest at `path`.
    ///
    /// # Errors
    ///
    /// [`DataError::Manifest`] for unreadable, unparsable or invalid content.
    pub fn load(path: &Path) -> Result<Self, DataError> {
        let fail = |message: String| DataError::Manifest {
            path: path.to_owned(),
            message,
        };
        let file: File = read_document(path).map_err(|error| fail(error.to_string()))?;
        let mut names = std::collections::BTreeSet::new();
        for entry in &file.entry {
            check(entry).map_err(|message| fail(format!("entry `{}`: {message}", entry.name)))?;
            if !names.insert(entry.name.as_str()) {
                return Err(fail(format!("entry `{}` is declared twice", entry.name)));
            }
        }
        Ok(Self {
            path: path.to_owned(),
            entries: file.entry,
        })
    }

    /// A manifest over `entries` that skips file validation, for tests that serve
    /// a local address.
    #[cfg(test)]
    pub(super) fn from_entries(path: PathBuf, entries: Vec<Entry>) -> Self {
        Self { path, entries }
    }

    /// The manifest file.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The entries in file order.
    #[must_use]
    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// The entry called `name`.
    ///
    /// # Errors
    ///
    /// [`DataError::NoEntry`] when there is none.
    pub fn entry(&self, name: &str) -> Result<&Entry, DataError> {
        self.entries
            .iter()
            .find(|entry| entry.name == name)
            .ok_or_else(|| DataError::NoEntry {
                path: self.path.clone(),
                name: name.to_owned(),
            })
    }

    /// The directory holding the manifest, where `append` paths start.
    #[must_use]
    pub fn base(&self) -> &Path {
        match self.path.parent() {
            Some(parent) => parent,
            None => Path::new(""),
        }
    }
}

impl Entry {
    /// The cache file name, `<name>.<extension>`.
    #[must_use]
    pub fn file_name(&self) -> String {
        format!("{}.{}", self.name, self.extension)
    }
}

fn check(entry: &Entry) -> Result<(), String> {
    let word = |key: &str, value: &str| {
        let fits = !value.is_empty()
            && value
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.'));
        if fits {
            Ok(())
        } else {
            Err(format!("`{key}` must be letters, digits, `-`, `_` or `.`"))
        }
    };
    word("name", &entry.name)?;
    word("extension", &entry.extension)?;
    for (key, value) in [("version", &entry.version), ("licence", &entry.licence)] {
        if value.trim().is_empty() || value.contains('\n') {
            return Err(format!("`{key}` must be one non-empty line"));
        }
    }
    if let Some(source) = &entry.source
        && (source.trim().is_empty() || source.contains('\n'))
    {
        return Err("`source` must be one non-empty line".to_owned());
    }
    if !entry.url.starts_with("https://") {
        return Err("`url` must start with https://".to_owned());
    }
    if entry.sha256.len() != 64 || !entry.sha256.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err("`sha256` must be 64 hexadecimal digits".to_owned());
    }
    Ok(())
}

#[cfg(test)]
#[path = "manifest_tests.rs"]
mod tests;
