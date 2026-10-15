//! Lazy `BinaryCIF` document container.

mod document;

pub use document::BinaryDocument;
pub(crate) use document::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
