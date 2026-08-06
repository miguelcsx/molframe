//! The chains of a structure.
//!
//! Both identifier namespaces are stored. The file's own normalised label is
//! always present; the depositor's label is optional, and absent is a different
//! statement from equal-to-the-label. Choosing between them at read time is what
//! makes a result impossible to trace back to either the file or the paper, so
//! nothing here chooses.

use crate::index::{ChainIndex, EntityIndex};
use crate::optional::OptionalSymbol;
use crate::symbol::SymbolId;
use std::ops::Range;

/// What kind of polymer a chain is, where it is one.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
#[non_exhaustive]
pub enum PolymerKind {
    /// Not a polymer.
    #[default]
    None,
    /// A polypeptide.
    Protein,
    /// Deoxyribonucleic acid.
    Dna,
    /// Ribonucleic acid.
    Rna,
    /// A chain containing both deoxyribo- and ribonucleotides.
    NucleicHybrid,
    /// A polysaccharide.
    Saccharide,
    /// A polymer of some other kind.
    Other,
}

impl PolymerKind {
    /// Returns true for any polymer.
    #[must_use]
    pub const fn is_polymer(self) -> bool {
        !matches!(self, Self::None)
    }

    /// Returns true for a nucleic acid of either kind, or a hybrid.
    #[must_use]
    pub const fn is_nucleic(self) -> bool {
        matches!(self, Self::Dna | Self::Rna | Self::NucleicHybrid)
    }
}

/// Everything a chain is created with.
#[derive(Clone, Copy, Debug)]
pub struct ChainRecord {
    /// The normalised chain label.
    pub label_asym_id: SymbolId,
    /// The depositor's chain label, if the file carried one.
    pub auth_asym_id: OptionalSymbol,
    /// The species this chain instantiates.
    pub entity: EntityIndex,
    /// What kind of polymer this chain is, if any.
    pub polymer_kind: PolymerKind,
}

/// The chain table.
#[derive(Clone, Debug, Default)]
pub struct ChainTable {
    first_residue: Vec<u32>,
    residue_count: Vec<u32>,
    label_asym_id: Vec<SymbolId>,
    auth_asym_id: Vec<OptionalSymbol>,
    entity: Vec<EntityIndex>,
    polymer_kind: Vec<PolymerKind>,
}

impl ChainTable {
    /// The number of chains.
    #[must_use]
    pub fn len(&self) -> usize {
        self.label_asym_id.len()
    }

    /// Returns true when the structure has no chains.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.label_asym_id.is_empty()
    }

    /// Appends a chain covering a range of residues.
    pub fn push(&mut self, record: ChainRecord, residues: Range<u32>) -> ChainIndex {
        let position = self.label_asym_id.len() as u32;
        self.first_residue.push(residues.start);
        self.residue_count
            .push(residues.end.saturating_sub(residues.start));
        self.label_asym_id.push(record.label_asym_id);
        self.auth_asym_id.push(record.auth_asym_id);
        self.entity.push(record.entity);
        self.polymer_kind.push(record.polymer_kind);
        ChainIndex::new(position)
    }

    /// The residues this chain contains.
    #[must_use]
    pub fn residues(&self, chain: ChainIndex) -> Option<Range<u32>> {
        let first = *self.first_residue.get(chain.as_usize())?;
        let count = *self.residue_count.get(chain.as_usize())?;
        Some(first..first.saturating_add(count))
    }

    /// The normalised chain label.
    #[must_use]
    pub fn label_asym_id(&self, chain: ChainIndex) -> Option<SymbolId> {
        self.label_asym_id.get(chain.as_usize()).copied()
    }

    /// The depositor's chain label, if the file carried one.
    #[must_use]
    pub fn auth_asym_id(&self, chain: ChainIndex) -> Option<SymbolId> {
        self.auth_asym_id
            .get(chain.as_usize())
            .and_then(|symbol| symbol.get())
    }

    /// The species this chain instantiates.
    #[must_use]
    pub fn entity(&self, chain: ChainIndex) -> Option<EntityIndex> {
        self.entity.get(chain.as_usize()).copied()
    }

    /// What kind of polymer this chain is.
    #[must_use]
    pub fn polymer_kind(&self, chain: ChainIndex) -> Option<PolymerKind> {
        self.polymer_kind.get(chain.as_usize()).copied()
    }

    /// Every chain that instantiates `entity`.
    ///
    /// This is the copies-of-the-same-molecule question, answered by reading a
    /// column rather than by comparing sequences.
    pub fn instances_of(&self, entity: EntityIndex) -> impl Iterator<Item = ChainIndex> + '_ {
        self.entity
            .iter()
            .enumerate()
            .filter(move |(_, chain_entity)| **chain_entity == entity)
            .map(|(position, _)| ChainIndex::new(position as u32))
    }

    /// Every chain position.
    pub fn iter(&self) -> impl Iterator<Item = ChainIndex> + '_ {
        (0..self.label_asym_id.len() as u32).map(ChainIndex::new)
    }
}
