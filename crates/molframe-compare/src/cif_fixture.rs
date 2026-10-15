//! Test-only reader for inline CIF fixtures.

use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

/// Parses `text` as mmCIF, panicking with the findings when the fixture is invalid.
pub(crate) fn read_cif(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}
