//! First-pass metadata projection and coordinate-layout classification.

use super::keep_lowering_category;
use crate::document::{CifValueRef, Document};
use crate::lexer::Quoting;
use crate::parser::{ColumnId, DocumentSink, ValueSink};
use molframe_core::span::ByteSpan;
use std::collections::HashSet;

pub(super) struct Projection {
    pub(super) metadata: Document,
    pub(super) layout: AtomLayout,
}

#[derive(Default)]
pub(super) struct AtomLayout {
    pub(super) has_atom_site: bool,
    atom_loops: usize,
    atom_rows: usize,
    models: Vec<i64>,
    supported: bool,
}

impl AtomLayout {
    pub(super) fn supports_direct(&self) -> bool {
        self.supported && self.atom_loops == 1 && self.atom_rows > 0
    }

    pub(super) fn model_count(&self) -> usize {
        self.models.len()
    }
}

pub(super) struct ProjectionSink {
    metadata: DocumentSink,
    layout: AtomLayout,
    blocks: usize,
    in_atom_loop: bool,
    values_in_row: usize,
    row_model: Option<i64>,
    keep_category: fn(&str) -> bool,
    model_column: Option<ColumnId>,
}

impl ProjectionSink {
    pub(super) fn new(keep_category: fn(&str) -> bool) -> Self {
        Self {
            metadata: DocumentSink::with_filter(|_| true),
            layout: AtomLayout {
                supported: true,
                ..AtomLayout::default()
            },
            blocks: 0,
            in_atom_loop: false,
            values_in_row: 0,
            row_model: None,
            keep_category,
            model_column: None,
        }
    }

    fn in_first_block(&self) -> bool {
        self.blocks == 1
    }

    pub(super) const fn in_atom_loop(&self) -> bool {
        self.in_atom_loop
    }

    pub(super) fn metadata(&self) -> &Document {
        self.metadata.document()
    }
}

impl<'input> ValueSink<'input> for ProjectionSink {
    type Output = Projection;

    fn block(&mut self, name: &str) {
        self.blocks += 1;
        if self.in_first_block() {
            self.metadata.block(name);
        }
    }

    fn begin_loop(&mut self, tags: &[(Box<str>, Box<str>)]) {
        self.in_atom_loop = false;
        self.values_in_row = 0;
        self.row_model = None;
        self.model_column = None;
        if !self.in_first_block() || !tags.iter().any(|(category, _)| &**category == "atom_site") {
            return;
        }

        self.layout.has_atom_site = true;
        self.layout.atom_loops += 1;
        self.in_atom_loop = true;
        let mut items = HashSet::with_capacity(tags.len());
        for (category, item) in tags {
            if &**category == "atom_site" && !items.insert(&**item) {
                self.layout.supported = false;
            }
        }
        self.model_column = tags
            .iter()
            .position(|(category, item)| {
                category.as_ref() == "atom_site" && item.as_ref() == "pdbx_PDB_model_num"
            })
            .and_then(ColumnId::from_position);
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
        if !self.in_first_block() {
            return;
        }
        if keep_lowering_category(category) || (self.keep_category)(category) {
            self.metadata.value(category, item, text, quoting, span);
        }
        if category != "atom_site" {
            if self.in_atom_loop {
                self.values_in_row += 1;
            }
            return;
        }

        self.layout.has_atom_site = true;
        if !self.in_atom_loop {
            self.layout.supported = false;
            return;
        }
        self.values_in_row += 1;
        if self.model_column == Some(column) {
            self.row_model = CifValueRef::parse(text, quoting).as_integer();
        }
    }

    fn end_row(&mut self) {
        if !self.in_atom_loop || !self.in_first_block() {
            return;
        }
        self.layout.atom_rows += 1;
        let model = match self.row_model {
            Some(model) => model,
            None => 1,
        };
        if self.layout.models.last().copied() != Some(model) {
            self.layout.models.push(model);
        }
        self.values_in_row = 0;
        self.row_model = None;
    }

    fn end_loop(&mut self) {
        if self.in_atom_loop && self.values_in_row != 0 {
            self.layout.supported = false;
        }
        self.in_atom_loop = false;
        self.values_in_row = 0;
        self.row_model = None;
    }

    fn finish(self) -> Self::Output {
        Projection {
            metadata: self.metadata.finish(),
            layout: self.layout,
        }
    }
}
