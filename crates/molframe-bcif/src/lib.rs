//! `BinaryCIF` codecs, lazy document access, reading and writing.

#![forbid(unsafe_code)]

mod codec;
mod container;
mod encode;
mod messagepack;
mod reader;
mod writer;

pub use codec::{
    DataType, Decoded, DecodedStringColumn, DecodedStringIter, EncodedData, Encoding, decode,
    decode_f32_into,
};
pub use container::BinaryDocument;
pub use encode::{encode_floats, encode_integers, encode_interval, encode_strings};
pub use reader::{BcifBatchSource, BcifReader, read, read_document, read_with_document};
#[doc(hidden)]
pub use reader::{ProjectedReadResult, read_with_metadata, read_with_projection};
pub use writer::{
    write_document, write_structure, write_structure_to, write_structure_to_with_memory_limit,
    write_structure_with_options,
};
