//! Typed `ModelCIF` metadata and quality metrics over the shared mmCIF document.

#![forbid(unsafe_code)]

mod lower;
mod model;
mod write;

pub use lower::lower;
pub use model::*;
pub use write::{write_canonical, write_canonical_with_options};
