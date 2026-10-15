//! Lazy `BinaryCIF` document container.

mod document;
mod support;

pub use document::BinaryDocument;
pub(crate) use document::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
pub(crate) use support::{container_error, masked_value};
