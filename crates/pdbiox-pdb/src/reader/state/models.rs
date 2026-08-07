//! Multi-model state and dense-topology validation.
//!
//! Coordinate frames share storage only when atom identity and order agree.
//! The first model records a compact signature and every later row is compared
//! before it can enter a dense coordinate store.

use super::{ReadState, fixed};
use crate::reader::lines::Line;
use pdbiox_core::chunk::ChunkBuilder;
use pdbiox_core::diagnostic::{Code, Diagnostic, Severity};
use pdbiox_core::element::Element;
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::span::ByteSpan;
use pdbiox_core::symbol::{AltId, SymbolId};

/// Which part of a multi-model file is currently being consumed.
pub(super) enum ModelRead {
    Before,
    Reading,
    Ignoring,
}

/// Atom identity without coordinates or a model number.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) struct AtomSignature {
    pub(super) chain: SymbolId,
    pub(super) seq: i64,
    pub(super) ins_code: OptionalSymbol,
    pub(super) atom_name: SymbolId,
    pub(super) component_id: SymbolId,
    pub(super) alt_id: AltId,
    pub(super) element: Element,
}

impl ReadState<'_> {
    pub(super) fn model(&mut self, line: &Line<'_>) {
        let number =
            match fixed::integer(line.text, 11, 14).and_then(|number| i32::try_from(number).ok()) {
                Some(number) => number,
                None => self.model_number,
            };
        match self.model_read {
            ModelRead::Before => {}
            ModelRead::Reading if self.options.only_first_model => {
                self.model_read = ModelRead::Ignoring;
                return;
            }
            ModelRead::Reading => self.start_next_frame(),
            ModelRead::Ignoring => return,
        }
        self.model_number = number;
        self.model_numbers.push(number);
        self.model_read = ModelRead::Reading;
    }

    fn start_next_frame(&mut self) {
        self.close_chain();
        self.verify_frame_len();
        let builder = std::mem::replace(&mut self.builder, ChunkBuilder::new());
        let (chunks, coords) = builder.finish();
        if self.data.chunks.is_empty() {
            self.data.chunks = chunks.into();
        }
        self.frames.push(coords);
        self.topology_locked = true;
        self.residue_key = None;
        self.chain_label = None;
        self.atom_position = 0;
    }

    pub(super) fn observe_signature(&mut self, line: &Line<'_>, signature: AtomSignature) {
        if !self.topology_locked {
            self.signatures.push(signature);
            return;
        }
        if self.signatures.get(self.atom_position as usize) == Some(&signature) {
            return;
        }
        self.findings.push(
            Diagnostic::new(Code::E3010)
                .with_message("a later model does not describe the same atom topology")
                .with_severity(Severity::Breaking)
                .at(ByteSpan::empty(line.at)),
        );
    }

    pub(super) fn verify_frame_len(&mut self) {
        if self.topology_locked && self.atom_position as usize != self.signatures.len() {
            self.findings.push(
                Diagnostic::new(Code::E3010)
                    .with_message("a later model has a different atom count")
                    .with_severity(Severity::Breaking),
            );
        }
    }
}
