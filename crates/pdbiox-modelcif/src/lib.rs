//! Typed `ModelCIF` metadata and quality metrics over the shared mmCIF document.

#![forbid(unsafe_code)]

mod lower;
mod model;
mod read;
mod write;

pub use lower::lower;
pub use model::*;
pub use read::{ModelCifBatchSource, ModelCifProjection, read_compact, read_compact_with_options};
pub use write::{write_canonical, write_canonical_to, write_canonical_with_options};
