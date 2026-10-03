//! Strict loading of declarative JSON and TOML documents.
//!
//! One place decides how a file's syntax follows from its suffix, so a policy
//! file and a reference file cannot disagree about what `.TOML` means.

use serde::de::DeserializeOwned;
use std::path::Path;

/// Why a document could not be read into its schema.
#[derive(Debug, thiserror::Error)]
pub(crate) enum DocumentError {
    /// The file could not be read.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// JSON syntax or schema was invalid.
    #[error(transparent)]
    Json(#[from] serde_json::Error),
    /// TOML syntax or schema was invalid.
    #[error(transparent)]
    Toml(#[from] toml::de::Error),
    /// The suffix names neither JSON nor TOML.
    #[error("files must use a .json or .toml suffix")]
    UnsupportedFormat,
}

/// Reads `path` as JSON or TOML according to its suffix, case-insensitively.
pub(crate) fn read_document<T: DeserializeOwned>(path: &Path) -> Result<T, DocumentError> {
    let syntax = match path.extension().and_then(std::ffi::OsStr::to_str) {
        Some(extension) if extension.eq_ignore_ascii_case("json") => true,
        Some(extension) if extension.eq_ignore_ascii_case("toml") => false,
        _ => return Err(DocumentError::UnsupportedFormat),
    };
    let text = std::fs::read_to_string(path)?;
    if syntax {
        Ok(serde_json::from_str(&text)?)
    } else {
        Ok(toml::from_str(&text)?)
    }
}
