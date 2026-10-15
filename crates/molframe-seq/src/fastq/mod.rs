//! FASTQ sequence format.

mod codec;

pub use codec::{FastqError, FastqRecord, parse_fastq, write_fastq};
