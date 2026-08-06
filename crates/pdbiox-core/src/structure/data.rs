//! The normalised structure: what a file means, rather than how it was written.
//!
//! A structure is immutable and shared by reference. Editing produces a new one
//! rather than changing this one in place, which is what makes handing out a
//! pointer to the coordinates safe: a buffer given to a caller belongs to a
//! snapshot that can never change underneath it, and every derived index stays
//! valid for as long as the snapshot it was built against.
//!
//! Formatting is not preserved here. Reproducing whitespace, column widths and
//! quoting is the document's job; this layer preserves meaning.

use crate::chunk::AtomChunk;
use crate::coords::{CoordinateBlock, CoordinateGeneration};
use crate::index::ModelIndex;
use crate::symbol::Interner;
use crate::topology::Topology;
use std::sync::Arc;

/// What the entry as a whole says about itself.
///
/// The identifier is a variable-length string. Four characters was a property of
/// one file format, not of the archive, and extended identifiers are already
/// being issued.
#[derive(Clone, Debug, Default)]
pub struct EntryMetadata {
    /// The entry identifier, as deposited.
    pub id: Option<Box<str>>,
    /// The entry title.
    pub title: Option<Box<str>>,
    /// The experimental method.
    pub method: Option<Box<str>>,
    /// Resolution in ångström, where the method reports one.
    pub resolution: Option<f32>,
}

/// The crystallographic cell.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct UnitCell {
    /// Cell edge lengths in ångström.
    pub lengths: [f64; 3],
    /// Cell angles in degrees.
    pub angles: [f64; 3],
}

impl UnitCell {
    /// Returns true when the cell is the placeholder some files write when they
    /// have no cell at all.
    ///
    /// A unit cube with right angles is what a structure determined without
    /// crystallography carries, and treating it as a real cell would produce
    /// periodic images of a molecule that was never in a crystal.
    #[must_use]
    pub fn is_placeholder(&self) -> bool {
        self.lengths
            .iter()
            .all(|length| (length - 1.0).abs() < 1e-6)
            && self.angles.iter().all(|angle| (angle - 90.0).abs() < 1e-6)
    }
}

/// Where the positions of each model live.
#[derive(Clone, Debug)]
pub enum CoordinateStore {
    /// One model.
    Single(CoordinateBlock),
    /// Several models over identical topology, so only positions differ.
    ///
    /// This covers ensembles and trajectory frames, where annotations are shared
    /// and repeating them per frame would multiply the memory by the frame count.
    Dense {
        /// One block per model, each with the same atom count.
        frames: Vec<CoordinateBlock>,
    },
    /// Models whose atom content genuinely differs.
    ///
    /// Padding these into one rectangular block would fabricate atoms that were
    /// never modelled and destroy the difference between absent and unrecorded,
    /// so they are held as separate structures instead.
    Ragged {
        /// One structure per model.
        models: Vec<Structure>,
    },
}

impl CoordinateStore {
    /// The number of models.
    #[must_use]
    pub fn model_count(&self) -> usize {
        match self {
            Self::Single(_) => 1,
            Self::Dense { frames } => frames.len(),
            Self::Ragged { models } => models.len(),
        }
    }

    /// The positions of one model, where they are held directly.
    ///
    /// A ragged ensemble answers `None`: its models are structures in their own
    /// right and are reached through [`CoordinateStore::Ragged`].
    #[must_use]
    pub fn block(&self, model: ModelIndex) -> Option<&CoordinateBlock> {
        match self {
            Self::Single(block) if model.get() == 0 => Some(block),
            Self::Dense { frames } => frames.get(model.as_usize()),
            // A ragged ensemble's models are structures in their own right, and
            // a single-model store has no second model to offer.
            Self::Single(_) | Self::Ragged { .. } => None,
        }
    }

    /// Returns true when every model shares one topology.
    #[must_use]
    pub const fn is_dense(&self) -> bool {
        matches!(self, Self::Single(_) | Self::Dense { .. })
    }
}

/// Everything a structure holds.
#[derive(Clone, Debug)]
pub struct StructureData {
    /// What the entry says about itself.
    pub entry: EntryMetadata,
    /// The hierarchy.
    pub topology: Topology,
    /// The atoms, in chunks.
    pub chunks: Vec<AtomChunk>,
    /// The positions.
    pub coords: CoordinateStore,
    /// The identifiers this structure interned.
    pub dictionary: Interner,
    /// The crystallographic cell, where the file carried one.
    pub cell: Option<UnitCell>,
    /// How many times the positions have changed.
    pub generation: CoordinateGeneration,
}

/// An immutable structure, shared by reference.
///
/// Cloning one is a reference count, not a copy, so passing a structure to
/// another thread or holding several versions at once is cheap.
///
/// # Examples
///
/// ```
/// use pdbiox_core::structure::{Structure, StructureData};
///
/// let structure = Structure::new(StructureData::empty());
/// assert_eq!(structure.atom_count(), 0);
/// assert_eq!(structure.model_count(), 1);
/// ```
#[derive(Clone, Debug)]
pub struct Structure(Arc<StructureData>);

impl StructureData {
    /// A structure with nothing in it.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            entry: EntryMetadata::default(),
            topology: Topology::default(),
            chunks: Vec::new(),
            coords: CoordinateStore::Single(CoordinateBlock::new()),
            dictionary: Interner::new(),
            cell: None,
            generation: CoordinateGeneration::INITIAL,
        }
    }
}

impl Structure {
    /// Wraps structure data for sharing.
    #[must_use]
    pub fn new(data: StructureData) -> Self {
        Self(Arc::new(data))
    }

    /// The data behind the reference.
    #[must_use]
    pub fn data(&self) -> &StructureData {
        &self.0
    }

    /// The shared handle, for a view that must outlive this reference.
    #[must_use]
    pub fn shared(&self) -> Arc<StructureData> {
        Arc::clone(&self.0)
    }

    /// The number of atoms in one model.
    #[must_use]
    pub fn atom_count(&self) -> u32 {
        self.0.chunks.last().map_or(0, |chunk| chunk.atoms().end)
    }

    /// The number of models.
    #[must_use]
    pub fn model_count(&self) -> usize {
        self.0.coords.model_count()
    }

    /// The number of chains.
    #[must_use]
    pub fn chain_count(&self) -> usize {
        self.0.topology.chains.len()
    }

    /// The number of residues.
    #[must_use]
    pub fn residue_count(&self) -> usize {
        self.0.topology.residues.len()
    }

    /// The number of distinct chemical species.
    #[must_use]
    pub fn entity_count(&self) -> usize {
        self.0.topology.entities.len()
    }

    /// The positions of the first model, contiguously.
    ///
    /// This slice is the buffer itself, not a copy of it, which is what makes
    /// handing coordinates to a numeric array free.
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        match self.0.coords.block(ModelIndex::new(0)) {
            Some(block) => block.as_slice(),
            None => &[],
        }
    }

    /// The positions of one model.
    #[must_use]
    pub fn model_positions(&self, model: ModelIndex) -> Option<&[[f32; 3]]> {
        self.0.coords.block(model).map(CoordinateBlock::as_slice)
    }

    /// How many times the positions have changed.
    #[must_use]
    pub fn generation(&self) -> CoordinateGeneration {
        self.0.generation
    }

    /// The string an interned identifier names.
    #[must_use]
    pub fn resolve(&self, symbol: crate::symbol::SymbolId) -> Option<&str> {
        self.0.dictionary.resolve(symbol)
    }
}

impl From<StructureData> for Structure {
    fn from(data: StructureData) -> Self {
        Self::new(data)
    }
}

#[cfg(test)]
#[path = "data_tests.rs"]
mod tests;
