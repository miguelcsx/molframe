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
mod pdbml;
mod reader;
mod small_cif;
mod write;

pub use document::{Category, CifValue, Column, DataBlock, Document};
pub use lower::lower;
pub use parser::{ParseResult, Rows, parse, split_tag};
pub use pdbml::{PdbmlError, PdbmlReadError, parse_pdbml_document, read_pdbml, write_pdbml};
pub use reader::{CifReader, read, read_with_document};
pub use small_cif::{
    SmallCifAtom, SmallCifBond, SmallCifDialect, SmallCifError, SmallCifOptions, SmallCifStructure,
    lower_small_cif, lower_small_cif_with_options,
};
pub use write::{
    CifWriteError, CifWriteOptions, quote_text, render_value, write_canonical,
    write_canonical_with_options, write_preserving,
};
