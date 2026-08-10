//! Multi-model state and dense-topology validation.
//!
//! Coordinate frames share storage only when atom identity and order agree.
//! The first model records a compact signature and every later row is compared
//! before it can enter a dense coordinate store.

use super::{ReadState, fixed, hybrid36};
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
    pub(super) fn atom_signature(
        &mut self,
        line: &Line<'_>,
        raw_name: &str,
        element: Element,
    ) -> Option<AtomSignature> {
        let chain = self.intern(fixed::text(line.text, 22, 22));
        let Some(seq) = hybrid36::decode(fixed::raw(line.text, 23, 26), 4) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("residue number could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            return None;
        };
        let ins_text = fixed::text(line.text, 27, 27);
        let ins_code = if ins_text.is_empty() {
            OptionalSymbol::NONE
        } else {
            OptionalSymbol::some(self.intern(ins_text))
        };
        let atom_name = self.intern(raw_name.trim());
        let component_id = self.intern(fixed::text(line.text, 18, 20));
        let alt = fixed::text(line.text, 17, 17);
        let alt_id = if alt.is_empty() {
            Some(AltId::BLANK)
        } else {
            AltId::labelled(self.intern(alt))
        };
        let Some(alt_id) = alt_id else {
            self.findings.push(
                Diagnostic::new(Code::E1901)
                    .with_message("alternate-location identifier exceeds its encoding"),
            );
            return None;
        };
        Some(AtomSignature {
            chain,
            seq,
            ins_code,
            atom_name,
            component_id,
            alt_id,
            element,
        })
    }

    pub(super) fn model(&mut self, line: &Line<'_>) {
        let mut number = 0;
        if let Some(parsed) =
            fixed::integer(line.text, 11, 14).and_then(|value| i32::try_from(value).ok())
        {
            number = parsed;
        } else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("MODEL serial is absent or outside the supported range")
                    .at(ByteSpan::empty(line.at)),
            );
        }
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
