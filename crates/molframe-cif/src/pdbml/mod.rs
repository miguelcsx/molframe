//! PDBML/XML conversion.

mod codec;

pub use codec::{PdbmlError, PdbmlReadError, parse_pdbml_document, read_pdbml, write_pdbml};
