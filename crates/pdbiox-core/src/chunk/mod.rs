//! Atoms in chunks: columnar within a chunk, chunked across the structure.
//!
//! A chunk is at once the unit of storage, of cache residency, of parallel work
//! and of summary statistics. Sized so its hot columns stay in second-level
//! cache while a kernel runs, it bounds the working set by the chunk rather than
//! by the structure — a distance calculation over a ribosome touches the same
//! amount of memory at a time as one over a lysozyme.

mod atom;
mod builder;
mod parent;
mod stats;

pub use atom::*;
pub use builder::*;
pub use parent::*;
pub use stats::*;
