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
use std::sync::Arc;

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
    first_residue: Arc<Vec<u32>>,
    residue_count: Arc<Vec<u32>>,
    label_asym_id: Arc<Vec<SymbolId>>,
    auth_asym_id: Arc<Vec<OptionalSymbol>>,
    entity: Arc<Vec<EntityIndex>>,
    polymer_kind: Arc<Vec<PolymerKind>>,
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
        Arc::make_mut(&mut self.first_residue).push(residues.start);
        Arc::make_mut(&mut self.residue_count).push(residues.end.saturating_sub(residues.start));
        Arc::make_mut(&mut self.label_asym_id).push(record.label_asym_id);
        Arc::make_mut(&mut self.auth_asym_id).push(record.auth_asym_id);
        Arc::make_mut(&mut self.entity).push(record.entity);
        Arc::make_mut(&mut self.polymer_kind).push(record.polymer_kind);
        ChainIndex::new(position)
    }

    /// The residues this chain contains.
    #[must_use]
    pub fn residues(&self, chain: ChainIndex) -> Option<Range<u32>> {
        stored_range(&self.first_residue, &self.residue_count, chain.as_usize())
    }

    /// Finds the chain whose contiguous range contains `residue`.
    ///
    /// Chain ranges are ordered, so this is logarithmic in the chain count.
    #[must_use]
    pub fn containing(&self, residue: u32) -> Option<ChainIndex> {
        let candidate = self
            .first_residue
            .partition_point(|first| *first <= residue)
            .checked_sub(1)?;
        let chain = ChainIndex::new(candidate as u32);
        self.residues(chain)?.contains(&residue).then_some(chain)
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
            .copied()
            .and_then(OptionalSymbol::get)
    }

    /// Replaces both chain namespaces with one explicit label.
    ///
    /// Renaming both prevents a writer from silently preferring the old
    /// depositor label over the requested new normalised label.
    pub fn rename(&mut self, chain: ChainIndex, label: SymbolId) -> bool {
        let position = chain.as_usize();
        if position >= self.label_asym_id.len() || position >= self.auth_asym_id.len() {
            return false;
        }
        let Some(normalised) = Arc::make_mut(&mut self.label_asym_id).get_mut(position) else {
            return false;
        };
        *normalised = label;
        let Some(depositor) = Arc::make_mut(&mut self.auth_asym_id).get_mut(position) else {
            return false;
        };
        *depositor = OptionalSymbol::some(label);
        true
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

    /// Reclassifies a chain after component chemistry has been resolved.
    pub fn set_polymer_kind(&mut self, chain: ChainIndex, kind: PolymerKind) -> bool {
        let Some(value) = Arc::make_mut(&mut self.polymer_kind).get_mut(chain.as_usize()) else {
            return false;
        };
        *value = kind;
        true
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

fn stored_range(starts: &[u32], counts: &[u32], index: usize) -> Option<Range<u32>> {
    let start = *starts.get(index)?;
    let count = *counts.get(index)?;

    Some(start..start.saturating_add(count))
}
