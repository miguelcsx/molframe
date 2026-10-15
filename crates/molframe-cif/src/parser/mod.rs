//! CIF token-to-document parsing.

mod document;

pub use document::{CifEventSink, CifScalar, ParseResult, Rows, parse, parse_events, split_tag};
pub(crate) use document::{ColumnId, DocumentSink, ValueSink, parse_into};
