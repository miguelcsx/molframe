//! Direct one-pass text-CIF structure reading.
//!
//! Metadata projection, borrowed coordinate lowering and dense/ragged
//! classification share one lexer traversal. The lossless coordinate document
//! remains available only through the explicit `read_with_document` API.

mod extra;
mod projection;
mod reader;
mod row;
mod stream;

pub(super) use reader::{keep_lowering_category, read, read_with_metadata, read_with_projection};
