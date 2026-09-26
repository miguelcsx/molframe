//! Reading and writing PDBx/mmCIF.
//!
//! The format the archive actually distributes, and the one every other
//! structural format is now defined against. Reading it happens in two steps
//! that are deliberately separate: a document that preserves what the file said,
//! and a structure that says what it means. Keeping them apart is what lets a
//! file survive a round trip through categories this library has no
//! interpretation for.

#![forbid(unsafe_code)]

mod batch;
pub mod document;
pub mod lexer;
mod lower;
mod parser;
mod pdbml;
mod reader;
mod small_cif;
mod write;

pub use batch::MmcifBatchSource;
pub use document::{Category, CifValue, Column, DataBlock, Document};
pub use lower::lower;
#[doc(hidden)]
pub use lower::{
    AtomSiteRow, AtomSiteRowSink, Field, lower_atom_site_with, lower_ragged_atom_site_with,
    lower_single_atom_site_with,
};
#[doc(hidden)]
pub use parser::{CifEventSink, CifScalar, parse_events};
pub use parser::{ParseResult, Rows, parse, split_tag};
pub use pdbml::{PdbmlError, PdbmlReadError, parse_pdbml_document, read_pdbml, write_pdbml};
pub use reader::{CifReader, read, read_with_document};
#[doc(hidden)]
pub use reader::{ProjectedReadResult, read_with_metadata, read_with_projection};
pub use small_cif::{
    SmallCifAtom, SmallCifBond, SmallCifDialect, SmallCifError, SmallCifOptions, SmallCifStructure,
    lower_small_cif, lower_small_cif_with_options,
};
#[doc(hidden)]
pub use write::{
    CanonicalAtomRow, CanonicalProjection, CanonicalValue, canonical_projection,
    declared_polymer_types,
};
pub use write::{
    CifWriteError, CifWriteOptions, CifWriteToError, quote_text, render_value, write_canonical,
    write_canonical_to, write_canonical_with_options, write_preserving, write_preserving_to,
    write_quoted,
};
