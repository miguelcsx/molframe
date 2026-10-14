//! The heavy atoms of the reference residues and what the model has for them.
//!
//! Atoms are stored residue by residue in one flat table, so a pair of atoms is
//! a pair of integers and a residue is a range. Hydrogens are left out: models
//! rarely carry them and a reference that does would otherwise score atoms the
//! model cannot have. Building the table costs one pass over the atoms of the
//! residues involved plus a scan of each model residue per reference atom.

use crate::CompareError;
use molframe_core::Element;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::{ResidueRef, Structure};
use std::collections::BTreeMap;
use std::ops::Range;

/// A reference residue and the model residue that answers to it.
///
/// A reference residue without a counterpart still belongs to the comparison:
/// every distance involving its atoms is tested and none is preserved.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResidueCorrespondence {
    /// The reference residue.
    pub reference: ResidueIndex,
    /// The model residue, if the model has one.
    pub model: Option<ResidueIndex>,
}

impl From<(ResidueIndex, ResidueIndex)> for ResidueCorrespondence {
    fn from((reference, model): (ResidueIndex, ResidueIndex)) -> Self {
        Self {
            reference,
            model: Some(model),
        }
    }
}

pub(super) struct Table<'a> {
    pub names: Vec<&'a str>,
    pub reference: Vec<[f32; 3]>,
    pub model: Vec<Option<[f32; 3]>>,
    pub owner: Vec<u32>,
    pub spans: Vec<Range<u32>>,
    pub components: Vec<&'a str>,
    chain: Vec<u32>,
    ordinal: Vec<u32>,
}

impl<'a> Table<'a> {
    pub(super) fn build(
        reference: &'a Structure,
        model: &'a Structure,
        pairs: &[ResidueCorrespondence],
    ) -> Result<Self, CompareError> {
        let chain_of = chain_ordinals(reference);
        let mut table = Self {
            names: Vec::new(),
            reference: Vec::new(),
            model: Vec::new(),
            owner: Vec::new(),
            spans: Vec::with_capacity(pairs.len()),
            components: Vec::with_capacity(pairs.len()),
            chain: Vec::with_capacity(pairs.len()),
            ordinal: Vec::with_capacity(pairs.len()),
        };
        for (entry, pair) in pairs.iter().enumerate() {
            let residue = reference
                .data()
                .residue(pair.reference)
                .ok_or(CompareError::InvalidMapping)?;
            let counterpart = match pair.model {
                Some(index) => Some(
                    model
                        .data()
                        .residue(index)
                        .ok_or(CompareError::InvalidMapping)?,
                ),
                None => None,
            };
            let entry = u32::try_from(entry).map_err(|_| CompareError::InvalidMapping)?;
            let start = table.len()?;
            table.push_residue(entry, residue, counterpart);
            table.spans.push(start..table.len()?);
            table.components.push(match residue.name() {
                Some(name) => name,
                None => "",
            });
            let chain = match chain_of.get(&pair.reference.get()) {
                Some(chain) => *chain,
                None => u32::MAX,
            };
            table.chain.push(chain);
            table.ordinal.push(pair.reference.get());
        }
        if table
            .reference
            .iter()
            .flatten()
            .chain(table.model.iter().flatten().flatten())
            .any(|value| !value.is_finite())
        {
            return Err(CompareError::InvalidScoreInput);
        }
        Ok(table)
    }

    /// The identity assignment, one entry per reference atom.
    pub(super) fn identity(&self) -> Vec<u32> {
        (0..self.reference.len())
            .filter_map(|index| u32::try_from(index).ok())
            .collect()
    }

    fn len(&self) -> Result<u32, CompareError> {
        u32::try_from(self.reference.len()).map_err(|_| CompareError::InvalidMapping)
    }

    fn push_residue(
        &mut self,
        entry: u32,
        residue: ResidueRef<'a>,
        counterpart: Option<ResidueRef<'a>>,
    ) {
        let candidates: Vec<(&str, [f32; 3])> = match counterpart {
            Some(other) => heavy_atoms(other),
            None => Vec::new(),
        };
        for (name, position) in heavy_atoms(residue) {
            self.names.push(name);
            self.reference.push(position);
            self.owner.push(entry);
            let found = candidates
                .iter()
                .find(|(candidate, _)| *candidate == name)
                .map(|(_, position)| *position);
            self.model.push(found);
        }
    }

    /// Whether the atoms of two entries are tested against each other.
    pub(super) fn testable(&self, first: u32, second: u32, separation: u32) -> bool {
        if first == second {
            return false;
        }
        let (first, second) = (first as usize, second as usize);
        if self.chain[first] != self.chain[second] {
            return true;
        }
        self.ordinal[first].abs_diff(self.ordinal[second]) > separation
    }
}

/// The named heavy atoms of a residue, keeping the first of any repeated name.
fn heavy_atoms(residue: ResidueRef<'_>) -> Vec<(&str, [f32; 3])> {
    let mut atoms: Vec<(&str, [f32; 3])> = Vec::new();
    for atom in residue.atoms() {
        if atom.element().is_some_and(Element::is_hydrogen) {
            continue;
        }
        let (Some(name), Some(position)) = (atom.name(), atom.position()) else {
            continue;
        };
        if atoms.iter().all(|(seen, _)| *seen != name) {
            atoms.push((name, position));
        }
    }
    atoms
}

fn chain_ordinals(structure: &Structure) -> BTreeMap<u32, u32> {
    let mut ordinals = BTreeMap::new();
    for (chain, handle) in structure.data().chains().enumerate() {
        let chain = match u32::try_from(chain) {
            Ok(chain) => chain,
            Err(_) => u32::MAX,
        };
        for residue in handle.residues() {
            ordinals.insert(residue.index().get(), chain);
        }
    }
    ordinals
}
