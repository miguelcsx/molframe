//! Turning a document into a structure.
//!
//! The document says what the file contained; the structure says what it means.
//! Everything interpretive happens on this side of the line, and every decision
//! that the data did not force is reported.

mod atoms;
mod bonds;
mod diagnostics;
mod ensemble;
mod entry;
mod keys;
mod metadata;
mod ragged;
mod references;
mod stream;

pub use atoms::{AtomSiteRow, AtomSiteRowSink, Field};
pub use entry::{
    lower, lower_atom_site_with, lower_ragged_atom_site_with, lower_single_atom_site_with,
};
pub(crate) use stream::{
    StreamFrameParts, StreamModelBuilder, StreamModelParts, StreamedModels, finish_streamed,
    share_model_topology,
};
