//! Turning a document into a structure.
//!
//! The document says what the file contained; the structure says what it means.
//! Everything interpretive happens on this side of the line, and every decision
//! that the data did not force is reported.

mod atoms;
mod entry;
mod keys;

pub use entry::lower;
