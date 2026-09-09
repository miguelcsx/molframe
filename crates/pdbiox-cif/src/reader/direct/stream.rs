//! One-pass projection, coordinate lowering and ensemble classification.

use super::projection::{Projection, ProjectionSink};
use super::row::AtomRow;
use crate::document::CifValueRef;
use crate::lexer::Quoting;
use crate::lower::{
    AtomSiteRow, Field, StreamFrameParts, StreamModelBuilder, StreamModelParts, StreamedModels,
    share_model_topology,
};
use crate::parser::{ColumnId, ValueSink};
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::io::ReadOptions;
use pdbiox_core::span::ByteSpan;

pub(super) struct DirectOutput {
    pub(super) projection: Projection,
    pub(super) models: Option<StreamedModels>,
    pub(super) setup_errors: Vec<Diagnostic>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum AtomLoopState {
    Searching,
    Active,
    Finished,
}

pub(super) struct DirectSink<'options, 'input> {
    projection: ProjectionSink,
    options: &'options ReadOptions,
    row: AtomRow<'input>,
    current_number: Option<i64>,
    lowerer: Option<StreamModelBuilder<'options>>,
    models: ModelCollector,
    setup_errors: Vec<Diagnostic>,
    atom_loop: AtomLoopState,
    has_unrepresented_items: bool,
    accepting_rows: bool,
    first_capacity: usize,
    atom_columns: Vec<Option<Field>>,
}

impl<'options> DirectSink<'options, '_> {
    pub(super) fn new(
        options: &'options ReadOptions,
        keep_category: fn(&str) -> bool,
        input_bytes: usize,
    ) -> Self {
        Self {
            projection: ProjectionSink::new(keep_category),
            options,
            row: AtomRow::new(),
            current_number: None,
            lowerer: None,
            models: ModelCollector::new(),
            setup_errors: Vec::new(),
            atom_loop: AtomLoopState::Searching,
            has_unrepresented_items: false,
            accepting_rows: true,
            first_capacity: estimate_atom_rows(input_bytes),
            atom_columns: Vec::new(),
        }
    }

    fn begin_atom_loop(&mut self, tags: &[(Box<str>, Box<str>)]) {
        if self.atom_loop != AtomLoopState::Searching {
            return;
        }
        self.atom_loop = AtomLoopState::Active;
        self.atom_columns.clear();
        self.atom_columns
            .extend(tags.iter().map(|(category, item)| {
                if category.as_ref() == "atom_site" {
                    Field::from_item(item)
                } else {
                    None
                }
            }));
        self.has_unrepresented_items = tags
            .iter()
            .any(|(category, item)| category.as_ref() == "atom_site" && !AtomRow::supports(item));
        match StreamModelBuilder::new(
            self.projection.metadata(),
            self.options,
            self.first_capacity,
        ) {
            Ok(lowerer) => self.lowerer = Some(lowerer),
            Err(errors) => self.setup_errors.extend(errors),
        }
    }

    fn end_atom_row(&mut self) {
        let number = match self.row.integer(Field::ModelNum) {
            Some(number) => number,
            None => 1,
        };
        if self.current_number != Some(number) {
            self.close_model();
            if self.options.only_first_model && self.models.has_model() {
                self.accepting_rows = false;
            }
            if self.accepting_rows {
                self.current_number = Some(number);
                let reference = self.models.reference();
                if let Some(lowerer) = &mut self.lowerer {
                    lowerer.start(number, reference);
                }
            }
        }
        if self.accepting_rows
            && let Some(lowerer) = &mut self.lowerer
        {
            lowerer.feed(&self.row);
        }
        self.row.advance();
    }

    fn close_model(&mut self) {
        let closed = self.lowerer.as_mut().and_then(StreamModelBuilder::close);
        if let Some((parts, matches_reference)) = closed {
            self.models
                .push(parts, matches_reference && !self.has_unrepresented_items);
        }
        self.current_number = None;
    }
}

impl<'input> ValueSink<'input> for DirectSink<'_, 'input> {
    type Output = DirectOutput;

    fn block(&mut self, name: &str) {
        self.projection.block(name);
    }

    fn begin_loop(&mut self, tags: &[(Box<str>, Box<str>)]) {
        self.projection.begin_loop(tags);
        if self.projection.in_atom_loop() {
            self.begin_atom_loop(tags);
        }
    }

    fn value(
        &mut self,
        column: ColumnId,
        category: &str,
        item: &str,
        text: &'input str,
        quoting: Quoting,
        span: ByteSpan,
    ) {
        self.projection
            .value(column, category, item, text, quoting, span);
        if self.atom_loop != AtomLoopState::Active {
            return;
        }
        let Some(field) = self.atom_columns.get(column.position()).copied().flatten() else {
            return;
        };
        let value = CifValueRef::parse(text, quoting);
        self.row.set(field, value);
    }

    fn end_row(&mut self) {
        self.projection.end_row();
        if self.atom_loop == AtomLoopState::Active {
            self.end_atom_row();
        }
    }

    fn end_loop(&mut self) {
        if self.atom_loop == AtomLoopState::Active {
            self.close_model();
            self.atom_loop = AtomLoopState::Finished;
        }
        self.projection.end_loop();
    }

    fn finish(mut self) -> Self::Output {
        self.close_model();
        DirectOutput {
            projection: self.projection.finish(),
            models: self.models.finish(),
            setup_errors: self.setup_errors,
        }
    }
}

struct ModelCollector {
    first: Option<StreamModelParts>,
    additional: Vec<StreamFrameParts>,
    ragged: Option<Vec<StreamModelParts>>,
}

impl ModelCollector {
    const fn new() -> Self {
        Self {
            first: None,
            additional: Vec::new(),
            ragged: None,
        }
    }

    const fn has_model(&self) -> bool {
        self.first.is_some()
    }

    fn reference(&self) -> Option<&pdbiox_core::structure::StructureData> {
        self.ragged
            .is_none()
            .then(|| self.first.as_ref().map(|first| &first.data))
            .flatten()
    }

    fn push(&mut self, parts: StreamModelParts, matches_reference: bool) {
        if let Some(ragged) = &mut self.ragged {
            ragged.push(parts);
            return;
        }
        if self.first.is_none() {
            self.first = Some(parts);
            return;
        }
        if matches_reference {
            self.additional.push(StreamFrameParts {
                number: parts.number,
                findings: parts.findings,
                coordinates: parts.coordinates,
            });
            return;
        }
        self.begin_ragged(parts);
    }

    fn begin_ragged(&mut self, candidate: StreamModelParts) {
        let Some(first) = self.first.take() else {
            return;
        };
        let mut ragged = Vec::with_capacity(self.additional.len().saturating_add(2));
        for frame in self.additional.drain(..) {
            ragged.push(share_model_topology(&first, frame));
        }
        ragged.insert(0, first);
        ragged.push(candidate);
        self.ragged = Some(ragged);
    }

    fn finish(mut self) -> Option<StreamedModels> {
        if let Some(ragged) = self.ragged.take() {
            return Some(StreamedModels::Ragged(ragged));
        }
        self.first.map(|first| StreamedModels::Dense {
            first: Box::new(first),
            additional: self.additional,
        })
    }
}

fn estimate_atom_rows(input_bytes: usize) -> usize {
    const LARGE_INPUT: usize = 64 * 1024 * 1024;
    const APPROXIMATE_ROW_BYTES: usize = 90;
    if input_bytes < LARGE_INPUT {
        return 0;
    }
    input_bytes / APPROXIMATE_ROW_BYTES
}
