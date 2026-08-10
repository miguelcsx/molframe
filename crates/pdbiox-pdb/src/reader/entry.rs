//! The entry point: bytes in, structure and findings out.

use super::ensemble::{ModelSpec, ragged_models};
use super::lines::Lines;
use super::state::ReadState;
use crate::fixed;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{Format, InputBuffer, ReadOptions, ReadResult, Reader};
use pdbiox_core::structure::{CoordinateStore, Structure, StructureData};

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

    match ragged_models(text, options.only_first_model) {
        Ok(Some(models)) => return read_ragged(text, options, &models, variant),
        Ok(None) => {}
        Err(finding) => return Err(vec![finding]),
    }

    let mut state = ReadState::new(options, variant);
    for line in Lines::new(text) {
        let line = line.map_err(|finding| vec![finding])?;
        state.line(&line);
    }
    state.finish()
}

fn read_ragged(
    text: &str,
    options: &ReadOptions,
    specs: &[ModelSpec],
    variant: Format,
) -> ReadResult {
    let mut models = Vec::with_capacity(specs.len());
    let mut findings = Vec::new();
    for spec in specs {
        let (model, model_findings) = read_selected_model(text, options, spec.ordinal, variant)?;
        models.push(model);
        findings.extend(model_findings);
    }

    let mut data = StructureData::empty();
    if let Some(first) = models.first() {
        data.entry = first.data().entry.clone();
        data.cell = first.data().cell;
        data.extensions = first.data().extensions.clone();
    }
    for spec in specs {
        if data.topology.models.push(spec.number, 0..0).is_err() {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
    data.coords = CoordinateStore::Ragged { models };
    options.finish(Structure::new(data), findings)
}

fn read_selected_model(
    text: &str,
    options: &ReadOptions,
    target: usize,
    variant: Format,
) -> ReadResult {
    let mut state = ReadState::new(options, variant);
    let mut ordinal = 0;
    let mut current = None;
    for line in Lines::new(text) {
        let line = line.map_err(|finding| vec![finding])?;
        match fixed::record(line.text) {
            "MODEL" => {
                current = Some(ordinal);
                if ordinal == target {
                    state.line(&line);
                }
                ordinal += 1;
            }
            "ENDMDL" => {
                if current == Some(target) {
                    state.line(&line);
                }
                current = None;
            }
            "ATOM" | "HETATM" | "TER" if current == Some(target) => state.line(&line),
            record if crate::header::is_metadata_record(record) || record == "CONECT" => {
                state.line(&line);
            }
            _ => {}
        }
    }
    state.finish()
}

#[cfg(test)]
#[path = "entry_tests.rs"]
mod tests;
