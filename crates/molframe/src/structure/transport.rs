//! Explicit transport identities take O(chains + residues) work and memory.

use super::Structure;
use crate::{Code, Diagnostic};
use molframe_core::{
    ChainIndex, ResidueIndex,
    optional::{OptionalI32, OptionalSymbol},
    structure::Structure as CoreStructure,
    topology::{ChainRecord, ChainTable, ResidueRecord, ResidueTable},
};

fn failure(reason: &str) -> Diagnostic {
    Diagnostic::new(Code::E3001).with_message(reason.to_owned())
}

impl Structure {
    /// Generates unambiguous label identities for transport, retaining author identities.
    ///
    /// Deposited label identifiers are replaced explicitly. Atom rows, coordinates,
    /// insertion codes and author identifiers are retained without per-atom copies.
    ///
    /// # Errors
    /// Returns a diagnostic when the hierarchy or identifier capacity is invalid.
    pub fn with_transport_identifiers(&self) -> Result<Self, Diagnostic> {
        let source = self.engine().data();
        let mut data = source.clone();
        let mut chains = ChainTable::default();
        let mut residues = ResidueTable::default();
        for position in 0..self.chain_count() {
            let index = ChainIndex::new(
                u32::try_from(position).map_err(|error| failure(&error.to_string()))?,
            );
            let range = source
                .topology
                .chains
                .residues(index)
                .ok_or_else(|| failure("chain lacks a residue range"))?;
            let label = data
                .dictionary
                .intern(&format!("chain_{position}"))
                .map_err(|error| failure(&error.to_string()))?;
            chains
                .push(
                    ChainRecord {
                        label_asym_id: label,
                        auth_asym_id: match source.topology.chains.auth_asym_id(index) {
                            Some(value) => OptionalSymbol::some(value),
                            None => OptionalSymbol::NONE,
                        },
                        entity: source
                            .topology
                            .chains
                            .entity(index)
                            .ok_or_else(|| failure("chain lacks an entity"))?,
                        polymer_kind: source
                            .topology
                            .chains
                            .polymer_kind(index)
                            .ok_or_else(|| failure("chain lacks a polymer classification"))?,
                    },
                    range.clone(),
                )
                .map_err(|error| failure(&error.to_string()))?;
            for (offset, position) in range.enumerate() {
                let index = ResidueIndex::new(position);
                let label_seq =
                    i32::try_from(offset + 1).map_err(|error| failure(&error.to_string()))?;
                let optional = |value| match value {
                    Some(value) => OptionalSymbol::some(value),
                    None => OptionalSymbol::NONE,
                };
                residues
                    .push(
                        ResidueRecord {
                            label_comp_id: source
                                .topology
                                .residues
                                .label_comp_id(index)
                                .ok_or_else(|| failure("residue lacks a component"))?,
                            auth_comp_id: optional(source.topology.residues.auth_comp_id(index)),
                            label_seq_id: OptionalI32::some(label_seq),
                            auth_seq_id: match source.topology.residues.auth_seq_id(index) {
                                Some(value) => OptionalI32::some(value),
                                None => OptionalI32::NONE,
                            },
                            ins_code: optional(source.topology.residues.ins_code(index)),
                            het: source.topology.residues.is_het(index),
                        },
                        source
                            .topology
                            .residues
                            .atoms(index)
                            .ok_or_else(|| failure("residue lacks atoms"))?,
                    )
                    .map_err(|error| failure(&error.to_string()))?;
            }
        }
        data.topology.chains = chains;
        data.topology.residues = residues;
        Ok(CoreStructure::new(data).into())
    }
}
