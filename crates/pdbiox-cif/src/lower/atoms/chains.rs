//! Building chain rows and resolving their entity references.
//!
//! Chain identity is read once when a chain opens. The hot atom path therefore
//! compares only the interned normalised label, while depositor labels and
//! entity relations remain available without being repeated per atom.

use super::{AtomBuilder, AtomSiteRow, Field};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::index::EntityIndex;
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::symbol::SymbolId;
use pdbiox_core::topology::{ChainRecord, EntityKind, PolymerKind};

impl AtomBuilder<'_> {
    pub(super) fn open_chain(&mut self, rows: &dyn AtomSiteRow, label: u32) {
        self.close_chain();
        self.chain = Some(label);
        self.chain_first_residue = self.residue_position;
        self.chain_auth = match rows.identifier(Field::AuthAsymId) {
            Some(auth) => OptionalSymbol::some(self.intern(&auth)),
            None => OptionalSymbol::NONE,
        };
        self.chain_entity = self.entity_for(rows, SymbolId::from_raw(label));
    }

    pub(super) fn close_chain(&mut self) {
        if let Some(label) = self.chain
            && self.residue_position > self.chain_first_residue
        {
            let Some(entity) = self.chain_entity.or_else(|| self.unknown_entity()) else {
                self.findings.push(Diagnostic::new(Code::E3001));
                self.chain = None;
                return;
            };
            let polymer_kind = match self.data.topology.entities.kind(entity) {
                Some(EntityKind::Polymer) => PolymerKind::Other,
                _ => PolymerKind::None,
            };
            if self
                .data
                .topology
                .chains
                .push(
                    ChainRecord {
                        label_asym_id: SymbolId::from_raw(label),
                        auth_asym_id: self.chain_auth,
                        entity,
                        polymer_kind,
                    },
                    self.chain_first_residue..self.residue_position,
                )
                .is_err()
            {
                self.findings.push(Diagnostic::new(Code::E3001));
            }
        }
        self.chain = None;
        self.chain_auth = OptionalSymbol::NONE;
        self.chain_entity = None;
    }

    /// Resolves a chain's declared entity without guessing from its sequence.
    fn entity_for(&mut self, rows: &dyn AtomSiteRow, label: SymbolId) -> Option<EntityIndex> {
        if let Some(entity_id) = rows.identifier(Field::LabelEntityId) {
            let entity_id = self.intern(&entity_id);
            if let Some(entity) = self.data.topology.entities.find_by_id(entity_id) {
                return Some(entity);
            }
            return match self.data.topology.entities.push(
                entity_id,
                EntityKind::Unknown,
                OptionalSymbol::NONE,
                &[],
            ) {
                Ok(entity) => Some(entity),
                Err(error) => {
                    self.findings.push(
                        Diagnostic::new(Code::E3001).with_context("cause", error.to_string()),
                    );
                    None
                }
            };
        }
        if let Some(mapping) = self
            .asym_entities
            .iter()
            .find(|mapping| mapping.asym_id == label)
        {
            return Some(mapping.entity);
        }
        self.unknown_entity()
    }

    /// The single unknown entity used when the source declares no relation.
    fn unknown_entity(&mut self) -> Option<EntityIndex> {
        if let Some(entity) = self
            .data
            .topology
            .entities
            .iter()
            .find(|entity| self.data.topology.entities.id(*entity).is_none())
        {
            Some(entity)
        } else {
            match self.data.topology.entities.push_without_id(
                EntityKind::Unknown,
                OptionalSymbol::NONE,
                &[],
            ) {
                Ok(entity) => Some(entity),
                Err(error) => {
                    self.findings.push(
                        Diagnostic::new(Code::E3001).with_context("cause", error.to_string()),
                    );
                    None
                }
            }
        }
    }
}
