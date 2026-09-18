//! Dense-model identity and deposited model numbers.
//!
//! The generic DOM path can retain compact signatures while building several
//! frames. The direct path stores no first-model side buffer: later rows compare
//! against the already encoded topology and chunks through shared backing data.

use super::{AtomBuilder, AtomSiteRow};
use crate::lower::diagnostics::at_source_row;
use molframe_core::chunk::{AtomChunk, AtomRecord, ChunkBuilder};
use molframe_core::diagnostic::{Code, Diagnostic, Severity};
use molframe_core::element::Element;
use molframe_core::optional::{OptionalI32, OptionalSymbol};
use molframe_core::structure::StructureData;
use molframe_core::symbol::{AltId, SymbolId};
use molframe_core::topology::Topology;
use std::sync::Arc;

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

pub(super) struct ExpectedAtoms {
    topology: Topology,
    chunks: Arc<Vec<AtomChunk>>,
    chunk: usize,
    matches: bool,
}

impl ExpectedAtoms {
    fn new(reference: &StructureData) -> Self {
        Self {
            topology: reference.topology.clone(),
            chunks: Arc::clone(&reference.chunks),
            chunk: 0,
            matches: true,
        }
    }

    fn len(&self) -> usize {
        self.chunks
            .last()
            .map_or(0, |chunk| chunk.atoms().end as usize)
    }

    fn matches(&mut self, position: u32, signature: AtomSignature, record: &AtomRecord) -> bool {
        while self
            .chunks
            .get(self.chunk)
            .is_some_and(|chunk| position >= chunk.atoms().end)
        {
            self.chunk += 1;
        }
        let Some(chunk) = self.chunks.get(self.chunk) else {
            self.matches = false;
            return false;
        };
        let Some(local) = position.checked_sub(chunk.atoms().start) else {
            self.matches = false;
            return false;
        };
        let same = self.signature(chunk, local) == Some(signature)
            && same_record(chunk, local, &self.topology, record);
        self.matches &= same;
        same
    }

    fn signature(&self, chunk: &AtomChunk, local: u32) -> Option<AtomSignature> {
        let residue = chunk.residue(local, &self.topology.residues)?;
        let chain = self.topology.chains.containing(residue.get())?;
        Some(AtomSignature {
            chain: self.topology.chains.label_asym_id(chain)?.get(),
            label_seq: OptionalI32::from(self.topology.residues.label_seq_id(residue)),
            auth_seq: OptionalI32::from(self.topology.residues.auth_seq_id(residue)),
            ins_code: self
                .topology
                .residues
                .ins_code(residue)
                .map_or(crate::lower::keys::ResidueKey::ABSENT, SymbolId::get),
            atom_name: chunk.atom_name(local)?,
            auth_atom_name: OptionalSymbol::from(chunk.auth_atom_name(local)),
            component_id: chunk
                .alternate_component_id(local)
                .or_else(|| self.topology.residues.label_comp_id(residue))?,
            alt_id: chunk.alt_id(local)?,
            element: chunk.element(local)?,
        })
    }

    fn is_match(&self, observed: usize) -> bool {
        self.matches && observed == self.len()
    }
}

fn same_record(chunk: &AtomChunk, local: u32, topology: &Topology, candidate: &AtomRecord) -> bool {
    chunk.element(local) == Some(candidate.element)
        && chunk.atom_name(local) == Some(candidate.atom_name)
        && chunk.auth_atom_name(local) == candidate.auth_atom_name.get()
        && chunk.alternate_component_id(local) == candidate.alternate_component_id.get()
        && chunk.alt_id(local) == Some(candidate.alt_id)
        && chunk.residue(local, &topology.residues) == Some(candidate.residue)
        && chunk.occupancy(local) == Some(candidate.occupancy)
        && chunk.b_factor(local) == Some(candidate.b_factor)
        && chunk.formal_charge(local) == Some(candidate.formal_charge)
        && chunk.atom_site_id(local) == Some(candidate.atom_site_id)
}

impl AtomBuilder<'_> {
    /// Uses an externally validated single-model row stream.
    #[must_use]
    pub(crate) fn without_identity_tracking(mut self) -> Self {
        self.track_identity = false;
        self
    }

    /// Compares streamed rows against already materialised first-model columns.
    #[must_use]
    pub(crate) fn against(mut self, reference: &StructureData) -> Self {
        self.data.dictionary = reference.dictionary.clone();
        self.expected = Some(ExpectedAtoms::new(reference));
        self.track_identity = false;
        self
    }

    pub(crate) fn matches_expected(&self) -> bool {
        self.expected
            .as_ref()
            .is_none_or(|expected| expected.is_match(self.atom_position as usize))
    }

    /// Reserves the only atom-sized buffers whose extent is known by readers.
    pub(crate) fn reserve_atoms(&mut self, atoms: usize) {
        self.builder.reserve(atoms);
        if self.track_identity {
            self.signatures.reserve(atoms);
        }
    }

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

    pub(super) fn observe_atom(
        &mut self,
        rows: &dyn AtomSiteRow,
        signature: AtomSignature,
        record: &AtomRecord,
    ) {
        if let Some(expected) = &mut self.expected {
            let _ = expected.matches(self.atom_position, signature, record);
            return;
        }
        if !self.topology_locked {
            if self.track_identity {
                self.signatures.push(signature);
            }
            return;
        }
        if self.signatures.get(self.atom_position as usize) == Some(&signature) {
            return;
        }
        self.report_identity_mismatch(rows);
    }

    fn report_identity_mismatch(&mut self, rows: &dyn AtomSiteRow) {
        self.findings.push(at_source_row(
            Diagnostic::new(Code::E3010)
                .with_message("a later model does not describe the same atom topology")
                .with_severity(Severity::Breaking)
                .in_category("atom_site"),
            rows.row(),
        ));
    }

    pub(super) fn verify_frame_len(&mut self) {
        if self.expected.is_some() {
            return;
        }
        if self.topology_locked && self.atom_position as usize != self.signatures.len() {
            self.report_frame_len_mismatch();
        }
    }

    fn report_frame_len_mismatch(&mut self) {
        self.findings.push(
            Diagnostic::new(Code::E3010)
                .with_message("a later model has a different atom count")
                .with_severity(Severity::Breaking)
                .in_category("atom_site"),
        );
    }
}
