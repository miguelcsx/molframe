//! FASTA sequence format.

mod codec;

pub use codec::{FastaRecord, parse_fasta, write_fasta};
