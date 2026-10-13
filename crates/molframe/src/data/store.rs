//! Listing, verifying and fetching the entries of a manifest.
//!
//! Fetching is the only place the library touches the network, it happens only
//! on this explicit call, and it goes through the one verified download the
//! adapters crate provides: the bytes are accepted only when their SHA-256
//! matches the manifest. The result is converted, validated by the shared table
//! reader and written atomically into the cache directory, where every
//! computation finds it without further network access.
//!
//! Cost: one download bounded by the caller's byte ceiling, then linear passes
//! over the text; `list` and `verify` read local files only.

use super::convert::Converter;
use super::error::DataError;
use super::manifest::{Entry, Manifest};
use super::render::{bare_rows, read_back, table_text};
use molframe_adapters::{DownloadOptions, fetch_verified};
use std::path::{Path, PathBuf};

/// Where an entry stands in the cache.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Status {
    /// No cache file.
    Missing,
    /// The file declares the manifest's name, version and licence and was made
    /// from the manifest's pinned bytes.
    Current,
    /// The file is a valid table that no longer matches the manifest.
    Stale(String),
    /// The file is not a valid table.
    Invalid(String),
}

/// One entry and its cache state.
#[derive(Clone, Debug)]
pub struct Verified {
    /// The manifest entry.
    pub entry: Entry,
    /// Its cache file.
    pub path: PathBuf,
    /// Its state.
    pub status: Status,
}

/// A table that was just stored.
#[derive(Clone, Debug)]
pub struct Fetched {
    /// The cache file written.
    pub path: PathBuf,
    /// Bytes written.
    pub bytes: usize,
    /// Rows stored.
    pub rows: usize,
}

/// Every entry of `manifest` with the state of its cache file under `directory`.
///
/// Reads local files only.
#[must_use]
pub fn list(manifest: &Manifest, directory: &Path) -> Vec<Verified> {
    manifest
        .entries()
        .iter()
        .map(|entry| check(entry, directory))
        .collect()
}

/// Checks the cache files of `name`, or of every entry, against the manifest.
///
/// # Errors
///
/// [`DataError::NoEntry`] for a name the manifest lacks.
pub fn verify(
    manifest: &Manifest,
    directory: &Path,
    name: Option<&str>,
) -> Result<Vec<Verified>, DataError> {
    let chosen: Vec<&Entry> = match name {
        Some(name) => vec![manifest.entry(name)?],
        None => manifest.entries().iter().collect(),
    };
    Ok(chosen
        .into_iter()
        .map(|entry| check(entry, directory))
        .collect())
}

fn check(entry: &Entry, directory: &Path) -> Verified {
    let path = directory.join(entry.file_name());
    let status = status_of(entry, &path);
    Verified {
        entry: entry.clone(),
        path,
        status,
    }
}

fn status_of(entry: &Entry, path: &Path) -> Status {
    let text = match std::fs::read_to_string(path) {
        Ok(text) => text,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Status::Missing,
        Err(error) => return Status::Invalid(error.to_string()),
    };
    let stored = match read_back(&text) {
        Ok(stored) => stored,
        Err(error) => return Status::Invalid(error.to_string()),
    };
    if stored.rows == 0 {
        return Status::Invalid("the table has no rows".to_owned());
    }
    let wanted = entry.sha256.to_ascii_lowercase();
    let differs = [
        ("name", stored.name.as_str(), entry.name.as_str()),
        ("version", stored.version.as_str(), entry.version.as_str()),
        ("licence", stored.licence.as_str(), entry.licence.as_str()),
    ]
    .into_iter()
    .find(|(_, have, want)| have != want);
    if let Some((key, have, want)) = differs {
        return Status::Stale(format!("{key} is `{have}`, the manifest says `{want}`"));
    }
    match stored.digest {
        Some(digest) if digest != wanted => Status::Stale(
            "it was made from different source bytes than the manifest pins".to_owned(),
        ),
        _ => Status::Current,
    }
}

/// Downloads `name`, converts and validates it, and stores it under `directory`.
///
/// # Errors
///
/// A [`DataError`] for an unknown entry or converter, an unreadable append file,
/// a failed or mismatched download, a conversion the shared reader refuses, an
/// existing file when `replace` is false, or a write failure.
pub fn fetch(
    manifest: &Manifest,
    name: &str,
    directory: &Path,
    replace: bool,
    options: DownloadOptions,
) -> Result<Fetched, DataError> {
    let entry = manifest.entry(name)?;
    let convert = |message: String| DataError::Convert {
        name: entry.name.clone(),
        message,
    };
    let converter = match &entry.converter {
        Some(id) => Some(Converter::named(id).map_err(&convert)?),
        None => None,
    };
    let target = directory.join(entry.file_name());
    if target.exists() && !replace {
        return Err(DataError::Exists { path: target });
    }
    let mut appended = Vec::new();
    for relative in &entry.append {
        let path = manifest.base().join(relative);
        let text = std::fs::read_to_string(&path).map_err(|source| DataError::Read {
            path: path.clone(),
            source,
        })?;
        appended.extend(bare_rows(&text));
    }
    let download = fetch_verified(&entry.url, &entry.sha256, options).map_err(|source| {
        DataError::Download {
            name: entry.name.clone(),
            source,
        }
    })?;
    let text = String::from_utf8(download.bytes)
        .map_err(|error| convert(format!("the download is not UTF-8 text: {error}")))?;
    let stored = match converter {
        Some(converter) => {
            let mut rows = converter.rows(&text).map_err(&convert)?;
            rows.extend(appended);
            table_text(entry, &rows)
        }
        None if appended.is_empty() => text,
        None => return Err(convert("`append` needs a converter".to_owned())),
    };
    let check = read_back(&stored).map_err(|error| convert(error.to_string()))?;
    if check.rows == 0 || check.name != entry.name {
        return Err(convert(format!(
            "the stored table must declare the name `{}` and have rows",
            entry.name
        )));
    }
    write_atomic(&target, stored.as_bytes())?;
    Ok(Fetched {
        path: target,
        bytes: stored.len(),
        rows: check.rows,
    })
}

fn write_atomic(target: &Path, bytes: &[u8]) -> Result<(), DataError> {
    let fail = |source| DataError::Write {
        path: target.to_owned(),
        source,
    };
    if let Some(parent) = target
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).map_err(fail)?;
    }
    let mut partial = target.as_os_str().to_owned();
    partial.push(".part");
    let partial = PathBuf::from(partial);
    std::fs::write(&partial, bytes).map_err(fail)?;
    std::fs::rename(&partial, target).map_err(fail)
}

#[cfg(test)]
#[path = "store_tests.rs"]
mod tests;
