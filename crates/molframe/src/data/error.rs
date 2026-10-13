//! The one error family of manifest-driven table acquisition.
//!
//! Each variant maps to a registry code: a manifest or converted table that does
//! not validate is `E2201`, a file or download that cannot be read `E7101`, an
//! entry the manifest does not offer `E7102`, a cache that cannot be written
//! `E7901`.

use molframe_core::{Code, diagnostic_from};
use std::path::PathBuf;

/// Why a manifest operation failed.
#[derive(Debug, thiserror::Error)]
pub enum DataError {
    /// The manifest is unreadable or violates its schema.
    #[error("manifest {path}: {message}")]
    Manifest {
        /// The manifest file.
        path: PathBuf,
        /// What is wrong.
        message: String,
    },
    /// The manifest has no entry of that name.
    #[error("manifest {path} has no entry named `{name}`")]
    NoEntry {
        /// The manifest file.
        path: PathBuf,
        /// The requested name.
        name: String,
    },
    /// A file the entry needs could not be read.
    #[error("cannot read {path}: {source}")]
    Read {
        /// The file.
        path: PathBuf,
        /// The cause.
        source: std::io::Error,
    },
    /// The download failed or did not match its pinned digest.
    #[error("download of `{name}` failed: {source}")]
    Download {
        /// The entry.
        name: String,
        /// The cause.
        source: molframe_adapters::DownloadError,
    },
    /// The downloaded bytes did not convert into a valid table.
    #[error("entry `{name}`: {message}")]
    Convert {
        /// The entry.
        name: String,
        /// What is wrong.
        message: String,
    },
    /// The cache entry already exists and replacement was not requested.
    #[error("{path} already exists; replace it explicitly")]
    Exists {
        /// The cache file.
        path: PathBuf,
    },
    /// The cache file could not be written.
    #[error("cannot write {path}: {source}")]
    Write {
        /// The cache file.
        path: PathBuf,
        /// The cause.
        source: std::io::Error,
    },
}

diagnostic_from!(DataError, |error| match error {
    DataError::Manifest { .. } | DataError::Convert { .. } => Code::E2201,
    DataError::NoEntry { .. } => Code::E7102,
    DataError::Read { .. } | DataError::Download { .. } => Code::E7101,
    DataError::Exists { .. } | DataError::Write { .. } => Code::E7901,
});
