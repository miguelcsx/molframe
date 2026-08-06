//! The four tables together.

use super::chain::ChainTable;
use super::entity::EntityTable;
use super::model::ModelTable;
use super::residue::ResidueTable;
use crate::index::ResidueIndex;

/// The whole hierarchy.
#[derive(Clone, Debug, Default)]
pub struct Topology {
    /// Models, in deposition order.
    pub models: ModelTable,
    /// Chains, contiguous within their model.
    pub chains: ChainTable,
    /// Residues, contiguous within their chain.
    pub residues: ResidueTable,
    /// Chemical species, referenced by chains rather than containing them.
    pub entities: EntityTable,
}

impl Topology {
    /// The number of atoms the residue table accounts for.
    #[must_use]
    pub fn atom_count(&self) -> u32 {
        let Some(last_position) = self.residues.len().checked_sub(1) else {
            return 0;
        };

        let last = ResidueIndex::new(last_position as u32);

        match self.residues.atoms(last) {
            Some(range) => range.end,
            None => 0,
        }
    }
}

#[cfg(test)]
#[path = "hierarchy_tests.rs"]
mod tests;
