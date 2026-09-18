//! MMTF `MessagePack` schema, codecs, lowering and writing.

mod codec;
mod model;
mod reader;
mod schema;
mod writer;

pub use model::{
    MMTF_METADATA_EXTENSION, MmtfEntityMetadata, MmtfGroupMetadata, MmtfMetadata, MmtfOptionalField,
};
pub use reader::read_mmtf;
pub use writer::{write_mmtf, write_mmtf_to};

#[cfg(test)]
#[path = "writer_tests.rs"]
mod writer_tests;
