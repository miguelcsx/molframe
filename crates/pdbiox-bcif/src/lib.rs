//! `BinaryCIF` codecs, lazy document access, reading and writing.

#![forbid(unsafe_code)]

mod codec;
mod container;
mod encode;
mod reader;
mod writer;

pub use codec::{DataType, Decoded, EncodedData, Encoding, decode};
pub use container::BinaryDocument;
pub use encode::{encode_floats, encode_integers, encode_interval, encode_strings};
pub use reader::{BcifReader, read, read_document, read_with_document};
pub use writer::{write_document, write_structure, write_structure_with_options};
