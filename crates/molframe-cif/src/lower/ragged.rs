//! One-pass lowering of models with independent atom topology.

use super::atoms::{AtomBuilder, AtomSiteRow, AtomSiteRowSink, Field};
use super::entry::{finish_model, prepare_model};
use crate::document::DataBlock;
use molframe_core::diagnostic::Diagnostics;
use molframe_core::io::ReadOptions;
use molframe_core::structure::Structure;

pub(super) struct RaggedParts {
    pub(super) models: Vec<Structure>,
    pub(super) model_numbers: Vec<i64>,
    pub(super) findings: Diagnostics,
}

pub(super) struct RaggedBuilder<'a> {
    block: &'a DataBlock,
    options: &'a ReadOptions,
    current_number: Option<i64>,
    current: Option<AtomBuilder<'a>>,
    models: Vec<Structure>,
    model_numbers: Vec<i64>,
    findings: Diagnostics,
}

impl<'a> RaggedBuilder<'a> {
    pub(super) fn new(
        block: &'a DataBlock,
        options: &'a ReadOptions,
        findings: Diagnostics,
        model_capacity: usize,
    ) -> Self {
        Self {
            block,
            options,
            current_number: None,
            current: None,
            models: Vec::with_capacity(model_capacity),
            model_numbers: Vec::with_capacity(model_capacity),
            findings,
        }
    }

    fn start_model(&mut self, number: i64) {
        self.close_model();
        let (data, findings, asym_entities) = prepare_model(self.block, self.options, Vec::new());
        self.current = Some(AtomBuilder::new(
            data,
            findings,
            self.options,
            asym_entities,
        ));
        self.current_number = Some(number);
        self.model_numbers.push(number);
    }

    fn close_model(&mut self) {
        let Some(builder) = self.current.take() else {
            return;
        };
        let (data, findings, coords) = builder.finish();
        let (model, model_findings) = finish_model(self.block, data, findings, coords);
        self.findings.extend(model_findings);
        self.models.push(model);
    }

    pub(super) fn finish(mut self) -> RaggedParts {
        self.close_model();
        RaggedParts {
            models: self.models,
            model_numbers: self.model_numbers,
            findings: self.findings,
        }
    }

    pub(super) fn abort(mut self) -> Diagnostics {
        if let Some(builder) = self.current.take() {
            self.findings.extend(builder.abort().finish());
        }
        self.findings
    }
}

impl AtomSiteRowSink for RaggedBuilder<'_> {
    fn feed(&mut self, row: &dyn AtomSiteRow) {
        let number = match row.integer(Field::ModelNum) {
            Some(number) => number,
            None => 1,
        };
        if self.current_number != Some(number) {
            self.start_model(number);
        }
        if let Some(builder) = &mut self.current {
            builder.feed(row);
        }
    }
}
