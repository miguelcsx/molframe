//! Deterministic `BinaryCIF` writing.

mod atoms;
mod column;
mod connections;
mod document;
mod metadata;
mod structure;

pub use document::write_document;
pub use structure::{
    write_structure, write_structure_to, write_structure_to_with_memory_limit,
    write_structure_with_options,
};
