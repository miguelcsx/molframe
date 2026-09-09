//! Neutral, one-pass topology transfer and verified byte acquisition.

#![forbid(unsafe_code)]

mod download;
mod topology;
mod topology_import;

pub use download::{DownloadError, DownloadOptions, VerifiedDownload, fetch_verified};

pub use topology::{MISSING_STRING, TopologyBatch, TopologyBatchError};
pub use topology_import::TopologyImportError;
