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
use crate::column::Presence;
use crate::element::Element;
use crate::index::{AtomIndex, ChainIndex, EntityIndex, ModelIndex, ResidueIndex};
use crate::symbol::{AltId, SymbolId};
use crate::topology::PolymerKind;
use std::ops::Range;

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
    const fn new(data: &'a StructureData, index: ModelIndex) -> Self {
        Self { data, index }
    }

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
        let data = self.data;
        let range = range_or_empty(self.data.topology.models.chains(self.index));

        range.map(move |position| ChainRef::new(data, ChainIndex::new(position)))
    }

    /// The chain with this label, in either namespace.
    #[must_use]
    pub fn chain(self, label: &str) -> Option<ChainRef<'a>> {
        let wanted = self.data.dictionary.get(label)?;

        self.chains().find(|chain| chain.has_label(wanted))
    }
}

impl<'a> ChainRef<'a> {
    const fn new(data: &'a StructureData, index: ChainIndex) -> Self {
        Self { data, index }
    }

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
        let data = self.data;
        let range = range_or_empty(self.data.topology.chains.residues(self.index));

        range.map(move |position| ResidueRef::new(data, ResidueIndex::new(position)))
    }

    /// The residue with this number, in whichever namespace carries it.
    #[must_use]
    pub fn residue(self, number: i32) -> Option<ResidueRef<'a>> {
        self.residues().find(|residue| residue.has_number(number))
    }

    fn has_label(self, wanted: SymbolId) -> bool {
        self.label_asym_id() == Some(wanted) || self.auth_asym_id() == Some(wanted)
    }
}

impl<'a> ResidueRef<'a> {
    const fn new(data: &'a StructureData, index: ResidueIndex) -> Self {
        Self { data, index }
    }

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
        self.label_comp_id()
            .and_then(|symbol| self.data.dictionary.resolve(symbol))
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
        self.data
            .topology
            .residues
            .ins_code(self.index)
            .and_then(|symbol| self.data.dictionary.resolve(symbol))
    }

    /// Whether the file recorded this residue as a heterogen.
    #[must_use]
    pub fn is_het(self) -> bool {
        self.data.topology.residues.is_het(self.index)
    }

    /// The atoms of this residue.
    pub fn atoms(self) -> impl Iterator<Item = AtomRef<'a>> {
        let data = self.data;
        let range = range_or_empty(self.data.topology.residues.atoms(self.index));

        range.map(move |position| AtomRef::new(data, AtomIndex::new(position)))
    }

    /// The atom with this name.
    #[must_use]
    pub fn atom(self, name: &str) -> Option<AtomRef<'a>> {
        let wanted = self.data.dictionary.get(name)?;

        self.atoms().find(|atom| atom.has_name(wanted))
    }

    fn has_number(self, number: i32) -> bool {
        self.auth_seq_id() == Some(number) || self.label_seq_id() == Some(number)
    }
}

impl<'a> AtomRef<'a> {
    const fn new(data: &'a StructureData, index: AtomIndex) -> Self {
        Self { data, index }
    }

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
        let atoms = chunk.atoms();

        if !atoms.contains(&position) {
            return None;
        }

        Some((chunk, position - atoms.start))
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
        self.name_symbol()
            .and_then(|symbol| self.data.dictionary.resolve(symbol))
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

        chunk.b_factor(local).and_then(recorded_value)
    }

    /// The occupancy, where the file recorded one.
    #[must_use]
    pub fn occupancy(self) -> Option<f32> {
        let (chunk, local) = self.located()?;

        chunk.occupancy(local).and_then(recorded_value)
    }

    /// The position in the first model, where one was recorded.
    #[must_use]
    pub fn position(self) -> Option<[f32; 3]> {
        let (chunk, local) = self.located()?;

        if !chunk.has_position(local) {
            return None;
        }

        self.data
            .coords
            .block(ModelIndex::new(0))?
            .as_slice()
            .get(self.index.as_usize())
            .copied()
    }

    /// The residue this atom belongs to.
    #[must_use]
    pub fn residue(self) -> Option<ResidueRef<'a>> {
        let (chunk, local) = self.located()?;
        let index = chunk.residue(local, &self.data.topology.residues)?;

        Some(ResidueRef::new(self.data, index))
    }

    fn has_name(self, wanted: SymbolId) -> bool {
        self.name_symbol() == Some(wanted)
    }
}

impl StructureData {
    /// A handle on one model.
    #[must_use]
    pub fn model(&self, model: ModelIndex) -> Option<ModelRef<'_>> {
        (model.as_usize() < self.topology.models.len()).then_some(ModelRef::new(self, model))
    }

    /// Every model.
    pub fn models(&self) -> impl Iterator<Item = ModelRef<'_>> {
        self.topology
            .models
            .iter()
            .map(move |index| ModelRef::new(self, index))
    }

    /// Every chain, across every model.
    pub fn chains(&self) -> impl Iterator<Item = ChainRef<'_>> {
        self.topology
            .chains
            .iter()
            .map(move |index| ChainRef::new(self, index))
    }

    /// Every residue, across every chain.
    pub fn residues(&self) -> impl Iterator<Item = ResidueRef<'_>> {
        let count = residue_count_as_u32(self.topology.residues.len());

        (0..count).map(move |position| ResidueRef::new(self, ResidueIndex::new(position)))
    }

    /// Every atom, in file order.
    pub fn atoms(&self) -> impl Iterator<Item = AtomRef<'_>> {
        (0..self.atom_count()).map(move |position| AtomRef::new(self, AtomIndex::new(position)))
    }
}

fn range_or_empty(range: Option<Range<u32>>) -> Range<u32> {
    match range {
        Some(range) => range,
        None => 0..0,
    }
}

fn recorded_value<T>(entry: (T, Presence)) -> Option<T> {
    let (value, presence) = entry;

    presence.is_present().then_some(value)
}

fn residue_count_as_u32(count: usize) -> u32 {
    match u32::try_from(count) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    }
}

#[cfg(test)]
#[path = "handle_tests.rs"]
mod tests;
