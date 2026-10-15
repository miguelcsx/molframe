//! Structure reading from CIF documents.

mod direct;
mod structure;

pub use structure::{
    CifReader, ProjectedReadResult, read, read_with_document, read_with_metadata,
    read_with_projection,
};
