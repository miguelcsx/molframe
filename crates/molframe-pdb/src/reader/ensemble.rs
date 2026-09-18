//! Building genuinely different models in one additional file pass.
//!
//! The normal reader discovers a topology difference while it is already
//! streaming the file. Ragged storage then needs one more pass over the atom
//! records, never one pass per model. Shared metadata and connectivity records
//! are retained separately because every independent model must receive them.

use super::lines::{Line, Lines};
use super::state::ReadState;
use crate::fixed;
use crate::header::is_metadata_record;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::{Format, ReadOptions, ReadResult};
use molframe_core::structure::{CoordinateStore, Structure, StructureData};

/// Non-coordinate records replayed into every independent model.
pub(super) struct CommonRecords<'a> {
    lines: Vec<Line<'a>>,
}

impl<'a> CommonRecords<'a> {
    pub(super) const fn new() -> Self {
        Self { lines: Vec::new() }
    }

    pub(super) fn observe(&mut self, line: &Line<'a>) {
        let record = fixed::record(line.text);
        if is_metadata_record(record) || record == "CONECT" {
            self.lines.push(*line);
        }
    }

    fn state<'options>(
        &self,
        options: &'options ReadOptions,
        variant: Format,
    ) -> ReadState<'options> {
        let mut state = ReadState::new(options, variant);
        for line in &self.lines {
            state.line(line);
        }
        state
    }
}

/// Materialises all ragged models during one traversal of their disjoint rows.
pub(super) fn read_ragged(
    text: &str,
    options: &ReadOptions,
    variant: Format,
    common: &CommonRecords<'_>,
) -> ReadResult {
    let mut models = Vec::new();
    let mut numbers = Vec::new();
    let mut findings = Vec::new();
    let mut current = None;

    for line in Lines::new(text) {
        let line = line.map_err(|finding| vec![finding])?;
        match fixed::record(line.text) {
            "MODEL" => {
                finish_model(&mut current, &mut models, &mut numbers, &mut findings)?;
                let mut state = common.state(options, variant);
                state.line(&line);
                current = Some(state);
            }
            "ENDMDL" => {
                if let Some(state) = current.as_mut() {
                    state.line(&line);
                }
                finish_model(&mut current, &mut models, &mut numbers, &mut findings)?;
            }
            "ATOM" | "HETATM" => {
                if current.is_none() {
                    current = Some(common.state(options, variant));
                }
                if let Some(state) = current.as_mut() {
                    state.line(&line);
                }
            }
            "TER" => {
                if let Some(state) = current.as_mut() {
                    state.line(&line);
                }
            }
            _ => {}
        }
    }
    finish_model(&mut current, &mut models, &mut numbers, &mut findings)?;
    finish_ensemble(models, numbers, findings, options)
}

fn finish_model(
    current: &mut Option<ReadState<'_>>,
    models: &mut Vec<Structure>,
    numbers: &mut Vec<i32>,
    findings: &mut Vec<Diagnostic>,
) -> Result<(), Vec<Diagnostic>> {
    let Some(state) = current.take() else {
        return Ok(());
    };
    let (model, model_findings) = state.finish()?;
    let number = match model
        .data()
        .models()
        .next()
        .and_then(molframe_core::structure::ModelRef::number)
    {
        Some(number) => number,
        None => 1,
    };
    models.push(model);
    numbers.push(number);
    findings.extend(model_findings);
    Ok(())
}

fn finish_ensemble(
    models: Vec<Structure>,
    numbers: Vec<i32>,
    mut findings: Vec<Diagnostic>,
    options: &ReadOptions,
) -> ReadResult {
    let mut data = StructureData::empty();
    if let Some(first) = models.first() {
        data.entry = first.data().entry.clone();
        data.cell = first.data().cell;
        data.extensions = first.data().extensions.clone();
    }
    for number in numbers {
        if data.topology.models.push(number, 0..0).is_err() {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
    data.coords = CoordinateStore::Ragged { models };
    options.finish(Structure::new(data), findings)
}

#[cfg(test)]
#[path = "ensemble_tests.rs"]
mod tests;
