//! Common declarative dispatch for sequence and tree formats.

mod codec;
mod types;

pub use codec::{SequenceReader, SequenceWriter, read_sequence, write_sequence};
pub use types::{SequenceDocument, SequenceFormat, SequenceFormatError};
