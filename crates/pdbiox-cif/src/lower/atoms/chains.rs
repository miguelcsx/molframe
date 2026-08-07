//! Building chain rows and resolving their entity references.
//!
//! Chain identity is read once when a chain opens. The hot atom path therefore
//! compares only the interned normalised label, while depositor labels and
//! entity relations remain available without being repeated per atom.

use super::AtomBuilder;
use crate::parser::Rows;
use pdbiox_core::index::EntityIndex;
use pdbiox_core::optional::OptionalSymbol;
use pdbiox_core::symbol::SymbolId;
use pdbiox_core::topology::{ChainRecord, EntityKind, PolymerKind};

impl AtomBuilder<'_> {
    pub(super) fn open_chain(&mut self, rows: &Rows<'_>, label: u32) {
        self.close_chain();
        self.chain = Some(label);
        self.chain_first_residue = self.residue_position;
        self.chain_auth = match rows.identifier("auth_asym_id") {
            Some(auth) => OptionalSymbol::some(self.intern(&auth)),
            None => OptionalSymbol::NONE,
        };
        self.chain_entity = Some(self.entity_for(rows, SymbolId::from_raw(label)));
    }

    pub(super) fn close_chain(&mut self) {
        if let Some(label) = self.chain
            && self.residue_position > self.chain_first_residue
        {
            let entity = match self.chain_entity {
                Some(entity) => entity,
                None => self.fallback_entity(),
            };
            let polymer_kind = match self.data.topology.entities.kind(entity) {
                Some(EntityKind::Polymer) => PolymerKind::Other,
                _ => PolymerKind::None,
            };
            self.data.topology.chains.push(
                ChainRecord {
                    label_asym_id: SymbolId::from_raw(label),
                    auth_asym_id: self.chain_auth,
                    entity,
                    polymer_kind,
                },
                self.chain_first_residue..self.residue_position,
            );
        }
        self.chain = None;
        self.chain_auth = OptionalSymbol::NONE;
        self.chain_entity = None;
    }

    /// Resolves a chain's declared entity without guessing from its sequence.
    fn entity_for(&mut self, rows: &Rows<'_>, label: SymbolId) -> EntityIndex {
        if let Some(entity_id) = rows.identifier("label_entity_id") {
            let entity_id = self.intern(&entity_id);
            if let Some(entity) = self.data.topology.entities.find_by_id(entity_id) {
                return entity;
            }
            return self.data.topology.entities.push(
                entity_id,
                EntityKind::Unknown,
                OptionalSymbol::NONE,
                &[],
            );
        }
        if let Some(mapping) = self
            .asym_entities
            .iter()
            .find(|mapping| mapping.asym_id == label)
        {
            return mapping.entity;
        }
        self.fallback_entity()
    }

    /// The single fallback species used only when the file declares no mapping.
    fn fallback_entity(&mut self) -> EntityIndex {
        let id = self.intern("1");
        match self.data.topology.entities.find_by_id(id) {
            Some(entity) => entity,
            None => {
                self.data
                    .topology
                    .entities
                    .push(id, EntityKind::Unknown, OptionalSymbol::NONE, &[])
            }
        }
    }
}
