//! The entry point: bytes in, structure and findings out.

use super::lines::Lines;
use super::state::ReadState;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{Format, InputBuffer, ReadOptions, ReadResult, Reader};

/// The reader for the legacy fixed-column format.
#[derive(Clone, Copy, Debug, Default)]
pub struct PdbReader;

impl Reader for PdbReader {
    const FORMAT: Format = Format::Pdb;

    fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
        read(input, options)
    }
}

/// Reads a structure from fixed-column text.
///
/// # Errors
///
/// Returns the findings that stopped the read: text that is not valid UTF-8, or
/// a file with no coordinate records in it at all.
pub fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    let Ok(text) = str::from_utf8(input.as_bytes()) else {
        return Err(vec![
            Diagnostic::new(Code::E1201).with_message("input is not valid text"),
        ]);
    };

    let mut state = ReadState::new(options);
    for line in Lines::new(text) {
        state.line(&line);
    }
    state.finish()
}

#[cfg(test)]
#[path = "entry_tests.rs"]
mod tests;
