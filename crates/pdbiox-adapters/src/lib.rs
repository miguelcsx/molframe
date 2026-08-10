//! Typed, one-pass projections for external structural-bioinformatics tools.

#![forbid(unsafe_code)]

mod download;
mod topology;

pub use download::{DownloadError, DownloadOptions, VerifiedDownload, fetch_verified};

pub use topology::{
    ExportAtom, ExportBond, ExportChain, ExportResidue, TopologyExport, TopologyExportError,
};
