//! `BinaryCIF` document and structure reading.

mod batch;
mod direct;
mod index;
mod structure;

pub use batch::BcifBatchSource;
pub use direct::{DecodedColumns, decode_columns};
pub use structure::{
    BcifReader, ProjectedReadResult, read, read_document, read_with_document, read_with_metadata,
    read_with_projection,
};
