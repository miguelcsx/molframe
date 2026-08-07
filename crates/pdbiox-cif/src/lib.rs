//! Reading and writing PDBx/mmCIF.
//!
//! The format the archive actually distributes, and the one every other
//! structural format is now defined against. Reading it happens in two steps
//! that are deliberately separate: a document that preserves what the file said,
//! and a structure that says what it means. Keeping them apart is what lets a
//! file survive a round trip through categories this library has no
//! interpretation for.

#![forbid(unsafe_code)]

pub mod document;
pub mod lexer;
mod lower;
mod parser;
mod reader;
mod write;
mod write_bonds;
mod write_references;

pub use document::{Category, CifValue, Column, DataBlock, Document};
pub use lower::lower;
pub use parser::{ParseResult, Rows, parse, split_tag};
pub use reader::{CifReader, read, read_with_document};
pub use write::{write_canonical, write_preserving};
