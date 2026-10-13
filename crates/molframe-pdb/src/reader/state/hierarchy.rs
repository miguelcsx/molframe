//! Chain and residue lifecycle for the fixed-column reader.
//!
//! A chain or residue is opened when an atom's identity differs from the one in
//! progress and closed, with its atom range finalised, the moment the next one
//! begins. Keeping the open/close pairs together is what guarantees that every
//! range written to the topology tables covers exactly the atoms read into it.

use super::super::lines::Line;
use super::rows::sequence_identifier;
use super::{ReadState, ResidueKey};
use crate::fixed;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::{EntityIndex, ResidueIndex};
use molframe_core::optional::OptionalSymbol;
use molframe_core::symbol::SymbolId;
use molframe_core::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};

impl ReadState<'_> {
    pub(super) fn begin_chain_if_new(&mut self, chain: SymbolId) {
        if self.chain_label == Some(chain) {
            return;
        }
        self.close_chain();
        self.chain_label = Some(chain);
        self.chain_first_residue = self.residue_position;
    }

    pub(super) fn begin_residue_if_new(
        &mut self,
        line: &Line<'_>,
        chain: SymbolId,
        seq: i64,
        ins_code: OptionalSymbol,
        component_id: SymbolId,
    ) {
        let key = ResidueKey {
            chain,
            seq,
            ins_code,
        };
        if self.residue_key.as_ref() == Some(&key) {
            return;
        }
        // The residue that was being filled ends where this one begins.
        self.close_residue();
        self.residue_key = Some(key);

        let het = fixed::record(line.text) == "HETATM";

        if self
            .data
            .topology
            .residues
            .push(
                ResidueRecord {
                    label_comp_id: component_id,
                    auth_comp_id: OptionalSymbol::some(component_id),
                    label_seq_id: sequence_identifier(seq),
                    auth_seq_id: sequence_identifier(seq),
                    ins_code,
                    het,
                },
                self.atom_position..self.atom_position,
            )
            .is_err()
        {
            self.findings.push(Diagnostic::new(Code::E3001));
            return;
        }
        self.residue_position += 1;
    }

    /// Closes the residue and chain being built, fixing their atom ranges.
    ///
    /// Once the topology is locked, later models only contribute coordinates.
    /// Their `TER` and `ENDMDL` records must leave the first model's ranges
    /// alone: `atom_position` then counts atoms of the current frame, not of
    /// the topology, and would rewrite a residue's extent backwards.
    pub(super) fn close_chain(&mut self) {
        if self.topology_locked {
            return;
        }
        self.close_residue();
        let Some(label_asym_id) = self.chain_label.take() else {
            return;
        };
        if self.residue_position == self.chain_first_residue {
            return;
        }
        let Ok(entity) = self.unclassified_entity() else {
            self.findings.push(Diagnostic::new(Code::E3001));
            return;
        };
        if self
            .data
            .topology
            .chains
            .push(
                ChainRecord {
                    label_asym_id,
                    auth_asym_id: OptionalSymbol::some(label_asym_id),
                    entity,
                    polymer_kind: PolymerKind::None,
                },
                self.chain_first_residue..self.residue_position,
            )
            .is_err()
        {
            self.findings.push(Diagnostic::new(Code::E3001));
        }
        self.residue_key = None;
    }

    /// Extends the residue being built to cover the atoms read into it.
    pub(super) fn close_residue(&mut self) {
        if self.residue_position == 0 {
            return;
        }
        let residue = ResidueIndex::new(self.residue_position - 1);
        let Some(range) = self.data.topology.residues.atoms(residue) else {
            return;
        };
        if let Err(error) = self
            .data
            .topology
            .residues
            .set_atoms(residue, range.start..self.atom_position)
        {
            self.findings
                .push(Diagnostic::new(Code::E3001).with_context("cause", error.to_string()));
        }
    }

    /// The single entity every chain of a fixed-column file shares.
    ///
    /// The format has no concept of a chemical species distinct from a chain, so
    /// there is nothing to distinguish here. Recording one entity rather than
    /// inventing several says exactly that.
    pub(super) fn unclassified_entity(&mut self) -> Result<EntityIndex, Diagnostic> {
        if let Some(entity) = self.entity {
            return Ok(entity);
        }
        let id = self.intern("1");
        let entity =
            self.data
                .topology
                .entities
                .push(id, EntityKind::Unknown, OptionalSymbol::NONE, &[]);
        let entity = entity.map_err(|error| {
            Diagnostic::new(Code::E3001).with_context("cause", error.to_string())
        })?;
        self.entity = Some(entity);
        Ok(entity)
    }

    pub(super) fn intern(&mut self, text: &str) -> SymbolId {
        if let Ok(symbol) = self.data.dictionary.intern(text) {
            symbol
        } else {
            self.findings.push(
                Diagnostic::new(Code::E1901).with_message("the identifier dictionary is full"),
            );
            SymbolId::from_raw(0)
        }
    }

    pub(super) fn component_text(&self, component: SymbolId) -> Box<str> {
        match self.data.dictionary.resolve(component) {
            Some(text) => text.into(),
            None => "".into(),
        }
    }
}
