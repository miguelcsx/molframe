//! Dense-model identity and deposited model numbers.
//!
//! A dense coordinate store is valid only when every frame describes the same
//! atoms in the same order. The first model records a compact integer signature;
//! later models are checked against it before their positions are accepted.

use super::AtomBuilder;
use crate::parser::Rows;
use pdbiox_core::chunk::ChunkBuilder;
use pdbiox_core::diagnostic::{Code, Diagnostic, Severity};
use pdbiox_core::element::Element;
use pdbiox_core::optional::{OptionalI32, OptionalSymbol};
use pdbiox_core::symbol::{AltId, SymbolId};

/// The semantic identity of one atom, excluding its coordinates and model.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) struct AtomSignature {
    pub(super) chain: u32,
    pub(super) label_seq: OptionalI32,
    pub(super) auth_seq: OptionalI32,
    pub(super) ins_code: u32,
    pub(super) atom_name: SymbolId,
    pub(super) auth_atom_name: OptionalSymbol,
    pub(super) component_id: SymbolId,
    pub(super) alt_id: AltId,
    pub(super) element: Element,
}

impl AtomBuilder<'_> {
    pub(super) fn start_model(&mut self, model: i64) {
        if self.model != i64::MIN {
            self.close_chain();
            self.close_residue();
            self.verify_frame_len();
            let builder = std::mem::replace(&mut self.builder, ChunkBuilder::new());
            let (chunks, coords) = builder.finish();
            if self.data.chunks.is_empty() {
                self.data.chunks = chunks.into();
            }
            self.frames.push(coords);
            self.topology_locked = true;
            self.atom_position = 0;
            self.current = None;
            self.chain = None;
        }
        self.model = model;
        let number = if let Ok(number) = i32::try_from(model) {
            number
        } else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("model number is outside the supported integer range")
                    .in_category("atom_site"),
            );
            i32::MAX
        };
        self.model_numbers.push(number);
    }

    pub(super) fn observe_signature(&mut self, rows: &Rows<'_>, signature: AtomSignature) {
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
                .in_category("atom_site")
                .at_row(rows.row() as u32),
        );
    }

    pub(super) fn verify_frame_len(&mut self) {
        if self.topology_locked && self.atom_position as usize != self.signatures.len() {
            self.findings.push(
                Diagnostic::new(Code::E3010)
                    .with_message("a later model has a different atom count")
                    .with_severity(Severity::Breaking)
                    .in_category("atom_site"),
            );
        }
    }
}
