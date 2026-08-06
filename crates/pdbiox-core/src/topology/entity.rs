//! The distinct chemical species a structure contains.
//!
//! An entity is not a level of the hierarchy. Several chains may be copies of
//! one, and no chain contains one, so it lives in its own table with a reference
//! from each chain. That is what makes "which chains are copies of the same
//! molecule?" a lookup rather than a sequence comparison, and what lets sequence
//! and source metadata survive a round trip instead of being repeated per chain
//! or dropped.
//!
//! The canonical sequence held here is the sequence that *should* be present.
//! The residues in a chain are what was actually modelled, and the difference
//! between them is the unmodelled region — which is data, not an absence.

use crate::index::EntityIndex;
use crate::optional::OptionalSymbol;
use crate::symbol::SymbolId;
use std::ops::Range;

/// What kind of chemical species an entity is.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
#[non_exhaustive]
pub enum EntityKind {
    /// A polymer: protein, nucleic acid, or a hybrid of the two.
    Polymer,
    /// A single chemical component that is not part of a polymer.
    NonPolymer,
    /// Water.
    Water,
    /// A branched species, such as an oligosaccharide.
    Branched,
    /// The file did not say.
    #[default]
    Unknown,
}

/// The entity table.
///
/// Canonical sequences share one pool rather than owning a collection each, so
/// a structure with three hundred entities holds one allocation for their
/// sequences instead of three hundred.
#[derive(Clone, Debug, Default)]
pub struct EntityTable {
    kind: Vec<EntityKind>,
    id: Vec<SymbolId>,
    description: Vec<OptionalSymbol>,
    sequence_span: Vec<Range<u32>>,
    sequence_pool: Vec<SymbolId>,
}

impl EntityTable {
    /// The number of entities.
    #[must_use]
    pub fn len(&self) -> usize {
        self.kind.len()
    }

    /// Returns true when the structure declares no entities.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.kind.is_empty()
    }

    /// Appends an entity, returning its position.
    pub fn push(
        &mut self,
        id: SymbolId,
        kind: EntityKind,
        description: OptionalSymbol,
        canonical_sequence: &[SymbolId],
    ) -> EntityIndex {
        let start = self.sequence_pool.len() as u32;

        self.sequence_pool.extend_from_slice(canonical_sequence);

        let end = self.sequence_pool.len() as u32;
        let position = self.kind.len() as u32;

        self.kind.push(kind);
        self.id.push(id);
        self.description.push(description);
        self.sequence_span.push(start..end);

        EntityIndex::new(position)
    }

    /// What kind of species this entity is.
    #[must_use]
    pub fn kind(&self, entity: EntityIndex) -> Option<EntityKind> {
        self.kind.get(entity.as_usize()).copied()
    }

    /// The identifier the file gave this entity.
    #[must_use]
    pub fn id(&self, entity: EntityIndex) -> Option<SymbolId> {
        self.id.get(entity.as_usize()).copied()
    }

    /// The entity's description, if the file carried one.
    #[must_use]
    pub fn description(&self, entity: EntityIndex) -> Option<SymbolId> {
        self.description
            .get(entity.as_usize())
            .copied()
            .and_then(OptionalSymbol::get)
    }

    /// The sequence this entity should have, one component per position.
    ///
    /// Compare it against the residues a chain actually holds to find what was
    /// not modelled.
    #[must_use]
    pub fn canonical_sequence(&self, entity: EntityIndex) -> &[SymbolId] {
        let Some(span) = self.sequence_span.get(entity.as_usize()) else {
            return &[];
        };

        sequence_slice(&self.sequence_pool, span)
    }

    /// Every entity position.
    pub fn iter(&self) -> impl Iterator<Item = EntityIndex> + '_ {
        (0..self.kind.len() as u32).map(EntityIndex::new)
    }
}

fn sequence_slice<'a>(pool: &'a [SymbolId], span: &Range<u32>) -> &'a [SymbolId] {
    let Ok(start) = usize::try_from(span.start) else {
        return &[];
    };
    let Ok(end) = usize::try_from(span.end) else {
        return &[];
    };

    match pool.get(start..end) {
        Some(sequence) => sequence,
        None => &[],
    }
}
