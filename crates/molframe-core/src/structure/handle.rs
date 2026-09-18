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
use super::handle_helpers::{range_or_empty, recorded_value};
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
    pub(crate) data: &'a StructureData,
    pub(crate) index: ChainIndex,
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

    /// The chain at one position within this model.
    #[must_use]
    pub fn chain_at(self, position: usize) -> Option<ChainRef<'a>> {
        self.chains().nth(position)
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

    /// The normalised chain label as text.
    #[must_use]
    pub fn label(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.label_asym_id()?)
    }

    /// The depositor's chain label, where the file carried one.
    #[must_use]
    pub fn auth_asym_id(self) -> Option<SymbolId> {
        self.data.topology.chains.auth_asym_id(self.index)
    }

    /// The depositor's chain label as text.
    #[must_use]
    pub fn auth_label(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.auth_asym_id()?)
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

    /// The residue at one position within this chain.
    #[must_use]
    pub fn residue_at(self, position: usize) -> Option<ResidueRef<'a>> {
        self.residues().nth(position)
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

    /// The depositor's component code, where it differs or was recorded.
    #[must_use]
    pub fn auth_comp_id(self) -> Option<SymbolId> {
        self.data.topology.residues.auth_comp_id(self.index)
    }

    /// The depositor's component code as text.
    #[must_use]
    pub fn auth_name(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.auth_comp_id()?)
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

    /// The atom at one position within this residue.
    #[must_use]
    pub fn atom_at(self, position: usize) -> Option<AtomRef<'a>> {
        self.atoms().nth(position)
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

    /// The depositor's atom name, where it differs or was recorded.
    #[must_use]
    pub fn auth_name_symbol(self) -> Option<SymbolId> {
        let (chunk, local) = self.located()?;
        chunk.auth_atom_name(local)
    }

    /// The depositor's atom name as text.
    #[must_use]
    pub fn auth_name(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.auth_name_symbol()?)
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

    /// The alternate-location label as text, absent for the shared blank form.
    #[must_use]
    pub fn alt_label(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.alt_id()?.symbol()?)
    }

    /// The component identity this atom carries.
    ///
    /// Most atoms inherit the residue's primary component. An alternate
    /// conformation may carry a different identity, which takes precedence.
    #[must_use]
    pub fn component_id(self) -> Option<SymbolId> {
        let (chunk, local) = self.located()?;
        if let Some(component) = chunk.alternate_component_id(local) {
            return Some(component);
        }
        self.residue()?.label_comp_id()
    }

    /// The atom's effective component identity as text.
    #[must_use]
    pub fn component_name(self) -> Option<&'a str> {
        self.data.dictionary.resolve(self.component_id()?)
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

    /// The formal charge, where the file recorded one.
    #[must_use]
    pub fn formal_charge(self) -> Option<i8> {
        let (chunk, local) = self.located()?;
        let (value, presence) = chunk.formal_charge(local)?;
        presence.is_present().then_some(value)
    }

    /// The file-local atom-site identifier.
    #[must_use]
    pub fn atom_site_id(self) -> Option<u32> {
        let (chunk, local) = self.located()?;
        chunk.atom_site_id(local)
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

    /// A handle on one chain.
    #[must_use]
    pub fn chain(&self, chain: ChainIndex) -> Option<ChainRef<'_>> {
        (chain.as_usize() < self.topology.chains.len()).then_some(ChainRef {
            data: self,
            index: chain,
        })
    }

    /// A handle on one residue.
    #[must_use]
    pub fn residue(&self, residue: ResidueIndex) -> Option<ResidueRef<'_>> {
        (residue.as_usize() < self.topology.residues.len()).then_some(ResidueRef {
            data: self,
            index: residue,
        })
    }

    /// Every chain, across every model.
    pub fn chains(&self) -> impl Iterator<Item = ChainRef<'_>> {
        self.topology
            .chains
            .iter()
            .map(move |index| ChainRef::new(self, index))
    }

    /// The chain with this label, in either namespace.
    #[must_use]
    pub fn chain_named(&self, label: &str) -> Option<ChainRef<'_>> {
        let wanted = self.dictionary.get(label)?;
        self.chains().find(|chain| {
            chain.label_asym_id() == Some(wanted) || chain.auth_asym_id() == Some(wanted)
        })
    }

    /// Every residue, across every chain.
    pub fn residues(&self) -> impl Iterator<Item = ResidueRef<'_>> {
        (0..self.topology.residues.len()).filter_map(move |position| {
            u32::try_from(position)
                .ok()
                .map(|position| ResidueRef::new(self, ResidueIndex::new(position)))
        })
    }

    /// Every atom, in file order.
    pub fn atoms(&self) -> impl Iterator<Item = AtomRef<'_>> {
        (0..self.atom_count()).map(move |position| AtomRef::new(self, AtomIndex::new(position)))
    }

    /// A handle on one atom.
    #[must_use]
    pub fn atom(&self, atom: AtomIndex) -> Option<AtomRef<'_>> {
        let count = self.chunks.last().map_or(0, |chunk| chunk.atoms().end);
        (atom.get() < count).then_some(AtomRef {
            data: self,
            index: atom,
        })
    }
}

#[cfg(test)]
#[path = "handle_tests.rs"]
mod tests;
