//! The entry point: bytes in, structure and findings out.

use super::ensemble::{CommonRecords, read_ragged};
use super::lines::Lines;
use super::state::ReadState;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::{Format, InputBuffer, ReadOptions, ReadResult, Reader};

/// The reader for the legacy fixed-column format.
#[derive(Clone, Copy, Debug, Default)]
pub struct PdbReader;

impl Reader for PdbReader {
    const FORMAT: Format = Format::Pdb;

    fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
        read_as(input, options, Format::Pdb)
    }
}

/// Reads a structure from fixed-column text.
///
/// # Errors
///
/// Returns the findings that stopped the read: text that is not valid UTF-8, or
/// a file with no coordinate records in it at all.
pub fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    read_as(input, options, Format::Pdb)
}

/// Reads PQR coordinates with partial charge and radius annotations.
///
/// # Errors
///
/// Returns findings when text, coordinate fields or variant values cannot be
/// interpreted under the requested parse mode.
pub fn read_pqr(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    read_as(input, options, Format::Pqr)
}

/// Reads PDBQT coordinates with partial charge and `AutoDock` type annotations.
///
/// # Errors
///
/// Returns findings when text, coordinate fields or variant values cannot be
/// interpreted under the requested parse mode.
pub fn read_pdbqt(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    read_as(input, options, Format::Pdbqt)
}

fn read_as(input: &InputBuffer, options: &ReadOptions, variant: Format) -> ReadResult {
    let Ok(text) = str::from_utf8(input.as_bytes()) else {
        return Err(vec![
            Diagnostic::new(Code::E1201).with_message("input is not valid text"),
        ]);
    };

    let mut common = CommonRecords::new();
    let mut state = ReadState::new(options, variant);
    for line in Lines::new(text) {
        let line = line.map_err(|finding| vec![finding])?;
        common.observe(&line);
        state.line(&line);
    }
    if state.requires_ragged() {
        drop(state);
        return read_ragged(text, options, variant, &common);
    }
    state.finish()
}

#[cfg(test)]
#[path = "entry_tests.rs"]
mod tests;
