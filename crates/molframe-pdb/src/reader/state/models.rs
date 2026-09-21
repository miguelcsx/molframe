//! Multi-model state and dense-topology validation.
//!
//! The first model builds the typed topology. Every later row is compared
//! directly with that materialised representation and contributes only one
//! coordinate to its frame. A mismatch selects a second, single-pass ragged
//! read; no per-atom signature is retained beside the topology.

use super::{ReadState, fixed, hybrid36};
use crate::reader::lines::Line;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::element::Element;
use molframe_core::optional::OptionalSymbol;
use molframe_core::span::ByteSpan;
use molframe_core::symbol::{AltId, SymbolId};

/// Which part of a multi-model file is currently being consumed.
pub(super) enum ModelRead {
    Before,
    Reading,
    Ignoring,
}

/// One parsed atom identity, held only while its row is processed.
#[derive(Clone, Copy)]
pub(super) struct AtomIdentity {
    pub(super) chain: SymbolId,
    pub(super) seq: i64,
    pub(super) ins_code: OptionalSymbol,
    pub(super) atom_name: SymbolId,
    pub(super) component_id: SymbolId,
    pub(super) alt_id: AltId,
    pub(super) element: Element,
    pub(super) het: bool,
}

impl ReadState<'_> {
    pub(super) fn atom_identity(
        &mut self,
        line: &Line<'_>,
        raw_name: &str,
        element: Element,
    ) -> Option<AtomIdentity> {
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
        Some(AtomIdentity {
            chain,
            seq,
            ins_code,
            atom_name,
            component_id,
            alt_id,
            element,
            het: fixed::record(line.text) == "HETATM",
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

    pub(in crate::reader) fn requires_ragged(&mut self) -> bool {
        self.verify_frame_len();
        self.requires_ragged
    }

    fn start_next_frame(&mut self) {
        self.close_chain();
        if self.topology_locked {
            self.finish_frame();
        } else {
            self.lock_topology();
        }
        self.residue_key = None;
        self.chain_label = None;
        self.atom_position = 0;
        self.expected_chunk = 0;
    }

    fn lock_topology(&mut self) {
        let builder = std::mem::take(&mut self.builder);
        let (chunks, coords) = builder.finish();
        let atom_capacity = match chunks.last() {
            Some(chunk) => chunk.atoms().end as usize,
            None => 0,
        };
        self.data.chunks = chunks.into();
        self.frames.push(coords);
        self.frame = CoordinateBlock::with_capacity(atom_capacity);
        self.topology_locked = true;
    }

    fn finish_frame(&mut self) {
        self.verify_frame_len();
        let atom_capacity = self.expected_atom_count() as usize;
        let next = CoordinateBlock::with_capacity(atom_capacity);
        self.frames.push(std::mem::replace(&mut self.frame, next));
    }

    pub(super) fn later_atom(&mut self, line: &Line<'_>, signature: AtomIdentity) {
        let position = self.position_of(line);
        let occupancy = self.optional_real(line, 55, 60, "occupancy", 1.0);
        let b_factor = self.optional_real(line, 61, 66, "B factor", 0.0);
        let serial = super::serial_of(line);
        let standard_matches = self.stored_atom_matches(signature, occupancy, b_factor, serial);
        let variant_matches = self.variant_atom_matches(line);
        self.requires_ragged |= !standard_matches || !variant_matches;

        self.frame.push(position_or_nan(position));
        self.atom_position += 1;
        self.saw_atoms = true;
    }

    fn stored_atom_matches(
        &mut self,
        signature: AtomIdentity,
        occupancy: (f32, molframe_core::column::Presence),
        b_factor: (f32, molframe_core::column::Presence),
        serial: u32,
    ) -> bool {
        while let Some(chunk) = self.data.chunks.get(self.expected_chunk) {
            if self.atom_position < chunk.atoms().end {
                break;
            }
            self.expected_chunk += 1;
        }
        let Some(chunk) = self.data.chunks.get(self.expected_chunk) else {
            return false;
        };
        let atoms = chunk.atoms();
        let Some(local) = self.atom_position.checked_sub(atoms.start) else {
            return false;
        };
        let Some(residue) = chunk.residue(local, &self.data.topology.residues) else {
            return false;
        };
        let Some(chain) = self.data.topology.chains.containing(residue.get()) else {
            return false;
        };
        let expected_component = match chunk.alternate_component_id(local) {
            Some(component) => Some(component),
            None => self.data.topology.residues.label_comp_id(residue),
        };
        let sequence = i32::try_from(signature.seq).ok();

        chunk.element(local) == Some(signature.element)
            && chunk.atom_name(local) == Some(signature.atom_name)
            && chunk.alt_id(local) == Some(signature.alt_id)
            && chunk.atom_site_id(local) == Some(serial)
            && chunk.occupancy(local) == Some(occupancy)
            && chunk.b_factor(local) == Some(b_factor)
            && expected_component == Some(signature.component_id)
            && self.data.topology.residues.auth_seq_id(residue) == sequence
            && self.data.topology.residues.ins_code(residue) == signature.ins_code.get()
            && self.data.topology.residues.is_het(residue) == signature.het
            && self.data.topology.chains.label_asym_id(chain) == Some(signature.chain)
    }

    pub(super) fn verify_frame_len(&mut self) {
        if self.topology_locked && self.atom_position != self.expected_atom_count() {
            self.requires_ragged = true;
        }
    }

    pub(super) fn finish_dense_frames(&mut self) {
        if self.topology_locked {
            self.finish_frame();
        }
    }

    fn expected_atom_count(&self) -> u32 {
        match self.data.chunks.last() {
            Some(chunk) => chunk.atoms().end,
            None => 0,
        }
    }
}

fn position_or_nan(position: Option<[f32; 3]>) -> [f32; 3] {
    let Some(position) = position else {
        return [f32::NAN; 3];
    };
    position
}
