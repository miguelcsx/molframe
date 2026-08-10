//! Reading and writing the legacy fixed-column structure format.
//!
//! The format predates the archive's current one and remains the most widely
//! read structural format there is. It is also the most constrained: every field
//! has a fixed width, and a structure large enough to overflow one cannot be
//! written honestly. This crate reads what the archive contains — including the
//! extended counting scheme that files past the field limits use — and refuses
//! to write anything it would have to truncate.

#![forbid(unsafe_code)]

pub mod fixed;
mod header;
pub mod hybrid36;
mod mmtf;
mod reader;
mod writer;

pub use header::{PDB_HEADERS_EXTENSION, PdbHeaderRecord, PdbHeaders, PdbHeadersExt};
pub use mmtf::{
    MMTF_METADATA_EXTENSION, MmtfEntityMetadata, MmtfGroupMetadata, MmtfMetadata,
    MmtfOptionalField, read_mmtf, write_mmtf,
};
pub use reader::{PdbReader, read, read_pdbqt, read_pqr};
pub use writer::{
    PdbIdentifierNamespace, PdbOptions, write, write_pdbqt, write_pqr, write_selected,
};
