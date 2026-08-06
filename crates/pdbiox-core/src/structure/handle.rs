//! Walking the hierarchy without building objects.
//!
//! A handle is a borrow of the structure plus a position — two words, `Copy`,
//! and allocation-free. Iterating a million atoms therefore allocates nothing at
//! all, which is the point on which an object-per-atom model becomes unusable at
//! scale while still reading naturally.
//!
//! The hierarchy is a computed view over the tables, not the storage. These
//! handles are how it is computed.

use super::data::StructureData;
use crate::chunk::AtomChunk;
use crate::element::Element;
use crate::index::{AtomIndex, ChainIndex, EntityIndex, ModelIndex, ResidueIndex};
use crate::symbol::{AltId, SymbolId};
use crate::topology::PolymerKind;

/// One model of a structure.
#[derive(Clone, Copy, Debug)]
pub struct ModelRef<'a> {
    data: &'a StructureData,
    index: ModelIndex,
}

/// One chain.
#[derive(Clone, Copy, Debug)]
pub struct ChainRef<'a> {
    data: &'a StructureData,
    index: ChainIndex,
}

/// One residue.
#[derive(Clone, Copy, Debug)]
pub struct ResidueRef<'a> {
    data: &'a StructureData,
    index: ResidueIndex,
}

/// One atom.
#[derive(Clone, Copy, Debug)]
pub struct AtomRef<'a> {
    data: &'a StructureData,
    index: AtomIndex,
}

impl<'a> ModelRef<'a> {
    /// This model's position.
    #[must_use]
    pub const fn index(self) -> ModelIndex {
        self.index
    }

    /// The number the model was deposited under, which is never reassigned.
    #[must_use]
    pub fn number(self) -> Option<i32> {
        self.data.topology.models.model_num(self.index)
    }

    /// The chains of this model.
    pub fn chains(self) -> impl Iterator<Item = ChainRef<'a>> {
        let range = match self.data.topology.models.chains(self.index) {
            Some(range) => range,
            None => 0..0,
        };
        let data = self.data;
        range.map(move |position| ChainRef {
            data,
            index: ChainIndex::new(position),
        })
    }

    /// The chain with this label, in either namespace.
    #[must_use]
    pub fn chain(self, label: &str) -> Option<ChainRef<'a>> {
        let wanted = self.data.dictionary.get(label)?;
        self.chains().find(|chain| {
            chain.label_asym_id() == Some(wanted) || chain.auth_asym_id() == Some(wanted)
        })
    }
}

impl<'a> ChainRef<'a> {
    /// This chain's position.
    #[must_use]
    pub const fn index(self) -> ChainIndex {
        self.index
    }

    /// The normalised chain label.
    #[must_use]
    pub fn label_asym_id(self) -> Option<SymbolId> {
        self.data.topology.chains.label_asym_id(self.index)
    }

    /// The depositor's chain label, where the file carried one.
    #[must_use]
    pub fn auth_asym_id(self) -> Option<SymbolId> {
        self.data.topology.chains.auth_asym_id(self.index)
    }

    /// The species this chain instantiates.
    #[must_use]
    pub fn entity(self) -> Option<EntityIndex> {
        self.data.topology.chains.entity(self.index)
    }

    /// What kind of polymer this chain is.
    #[must_use]
    pub fn polymer_kind(self) -> PolymerKind {
        match self.data.topology.chains.polymer_kind(self.index) {
            Some(kind) => kind,
            None => PolymerKind::None,
        }
    }

    /// The residues of this chain.
    pub fn residues(self) -> impl Iterator<Item = ResidueRef<'a>> {
        let range = match self.data.topology.chains.residues(self.index) {
            Some(range) => range,
            None => 0..0,
        };
        let data = self.data;
        range.map(move |position| ResidueRef {
            data,
            index: ResidueIndex::new(position),
        })
    }

    /// The residue with this number, in whichever namespace carries it.
    #[must_use]
    pub fn residue(self, number: i32) -> Option<ResidueRef<'a>> {
        self.residues().find(|residue| {
            residue.auth_seq_id() == Some(number) || residue.label_seq_id() == Some(number)
        })
    }
}

impl<'a> ResidueRef<'a> {
    /// This residue's position.
    #[must_use]
    pub const fn index(self) -> ResidueIndex {
        self.index
    }

    /// The normalised component code.
    #[must_use]
    pub fn label_comp_id(self) -> Option<SymbolId> {
        self.data.topology.residues.label_comp_id(self.index)
    }

    /// The component code as a string.
    #[must_use]
    pub fn name(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.label_comp_id()?)
    }

    /// The sequence position within the entity, absent for a non-polymer.
    #[must_use]
    pub fn label_seq_id(self) -> Option<i32> {
        self.data.topology.residues.label_seq_id(self.index)
    }

    /// The depositor's residue number.
    #[must_use]
    pub fn auth_seq_id(self) -> Option<i32> {
        self.data.topology.residues.auth_seq_id(self.index)
    }

    /// The insertion code, which is part of this residue's identity.
    #[must_use]
    pub fn ins_code(self) -> Option<&'a str> {
        let symbol = self.data.topology.residues.ins_code(self.index)?;
        self.data.dictionary.resolve(symbol)
    }

    /// Whether the file recorded this residue as a heterogen.
    #[must_use]
    pub fn is_het(self) -> bool {
        self.data.topology.residues.is_het(self.index)
    }

    /// The atoms of this residue.
    pub fn atoms(self) -> impl Iterator<Item = AtomRef<'a>> {
        let range = match self.data.topology.residues.atoms(self.index) {
            Some(range) => range,
            None => 0..0,
        };
        let data = self.data;
        range.map(move |position| AtomRef {
            data,
            index: AtomIndex::new(position),
        })
    }

    /// The atom with this name.
    #[must_use]
    pub fn atom(self, name: &str) -> Option<AtomRef<'a>> {
        let wanted = self.data.dictionary.get(name)?;
        self.atoms().find(|atom| atom.name_symbol() == Some(wanted))
    }
}

impl<'a> AtomRef<'a> {
    /// This atom's position in the flat atom order.
    #[must_use]
    pub const fn index(self) -> AtomIndex {
        self.index
    }

    /// The chunk holding this atom, and the position within it.
    ///
    /// Chunks tile the atom order in ascending ranges, so this is a binary
    /// search rather than a scan. It still runs once per attribute read, which
    /// is why a kernel walks chunks directly and leaves handles to the code
    /// that reads a few atoms rather than all of them.
    fn located(self) -> Option<(&'a AtomChunk, u32)> {
        let position = self.index.get();
        let candidate = self
            .data
            .chunks
            .partition_point(|chunk| chunk.atoms().end <= position);
        let chunk = self.data.chunks.get(candidate)?;
        chunk
            .atoms()
            .contains(&position)
            .then(|| (chunk, position - chunk.atoms().start))
    }

    /// The interned atom name.
    #[must_use]
    pub fn name_symbol(self) -> Option<SymbolId> {
        let (chunk, local) = self.located()?;
        chunk.atom_name(local)
    }

    /// The atom name as a string.
    #[must_use]
    pub fn name(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.name_symbol()?)
    }

    /// The element.
    #[must_use]
    pub fn element(self) -> Option<Element> {
        let (chunk, local) = self.located()?;
        chunk.element(local)
    }

    /// The alternate-location label.
    #[must_use]
    pub fn alt_id(self) -> Option<AltId> {
        let (chunk, local) = self.located()?;
        chunk.alt_id(local)
    }

    /// The temperature factor, where the file recorded one.
    #[must_use]
    pub fn b_factor(self) -> Option<f32> {
        let (chunk, local) = self.located()?;
        let (value, presence) = chunk.b_factor(local)?;
        presence.is_present().then_some(value)
    }

    /// The occupancy, where the file recorded one.
    #[must_use]
    pub fn occupancy(self) -> Option<f32> {
        let (chunk, local) = self.located()?;
        let (value, presence) = chunk.occupancy(local)?;
        presence.is_present().then_some(value)
    }

    /// The position in the first model, where one was recorded.
    #[must_use]
    pub fn position(self) -> Option<[f32; 3]> {
        let (chunk, local) = self.located()?;
        if !chunk.has_position(local) {
            return None;
        }
        let block = self.data.coords.block(ModelIndex::new(0))?;
        block.as_slice().get(self.index.as_usize()).copied()
    }

    /// The residue this atom belongs to.
    #[must_use]
    pub fn residue(self) -> Option<ResidueRef<'a>> {
        let (chunk, local) = self.located()?;
        let index = chunk.residue(local, &self.data.topology.residues)?;
        Some(ResidueRef {
            data: self.data,
            index,
        })
    }
}

impl StructureData {
    /// A handle on one model.
    #[must_use]
    pub fn model(&self, model: ModelIndex) -> Option<ModelRef<'_>> {
        (model.as_usize() < self.topology.models.len()).then_some(ModelRef {
            data: self,
            index: model,
        })
    }

    /// Every model.
    pub fn models(&self) -> impl Iterator<Item = ModelRef<'_>> {
        self.topology
            .models
            .iter()
            .map(move |index| ModelRef { data: self, index })
    }

    /// Every chain, across every model.
    pub fn chains(&self) -> impl Iterator<Item = ChainRef<'_>> {
        self.topology
            .chains
            .iter()
            .map(move |index| ChainRef { data: self, index })
    }

    /// Every residue, across every chain.
    pub fn residues(&self) -> impl Iterator<Item = ResidueRef<'_>> {
        (0..self.topology.residues.len() as u32).map(move |position| ResidueRef {
            data: self,
            index: ResidueIndex::new(position),
        })
    }

    /// Every atom, in file order.
    pub fn atoms(&self) -> impl Iterator<Item = AtomRef<'_>> {
        let count = self.chunks.last().map_or(0, |chunk| chunk.atoms().end);
        (0..count).map(move |position| AtomRef {
            data: self,
            index: AtomIndex::new(position),
        })
    }
}

#[cfg(test)]
#[path = "handle_tests.rs"]
mod tests;
