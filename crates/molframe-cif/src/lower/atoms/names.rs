//! The atom names already seen in the residue being filled.

use molframe_core::hashing::{IdentityBuildHasher, IdentityHashSet};
use molframe_core::symbol::{AltId, SymbolId};

/// The atom names already seen in the residue being filled.
///
/// A residue does not contain two atoms of the same name in the same alternate
/// location, so a repeat is what betrays a boundary the deposited identifiers
/// do not mark. Answering that question needs a set.
///
/// A linear scan is the right structure for an ordinary residue, which holds
/// about ten atoms — a set would cost more than it saves. It is the wrong
/// structure the moment a residue does not end: a coarse-grained model, a
/// lipid bilayer written as one residue, or any file whose annotations fail to
/// separate groups puts every atom in one residue, and the scan becomes
/// quadratic in the atom count. Both cases occur, so both are served: the scan
/// until it stops being cheap, an index after that.
#[derive(Default)]
pub(super) struct ResidueNames {
    names: Vec<(SymbolId, AltId)>,
    index: Option<IdentityHashSet<(u32, u32)>>,
}

impl ResidueNames {
    /// Atoms in a residue beyond which the linear scan stops paying.
    ///
    /// Well above any deposited residue, so an ordinary file never builds the
    /// index and never pays for it.
    const INDEX_THRESHOLD: usize = 64;

    pub(super) fn contains(&self, name: SymbolId, alt: AltId) -> bool {
        match &self.index {
            Some(index) => index.contains(&(name.get(), alt.get())),
            None => self.names.contains(&(name, alt)),
        }
    }

    pub(super) fn insert(&mut self, name: SymbolId, alt: AltId) {
        self.names.push((name, alt));
        match &mut self.index {
            Some(index) => {
                let _existing = index.insert((name.get(), alt.get()));
            }
            None if self.names.len() > Self::INDEX_THRESHOLD => {
                let mut index = IdentityHashSet::with_capacity_and_hasher(
                    self.names.len(),
                    IdentityBuildHasher::default(),
                );
                for &(name, alt) in &self.names {
                    let _existing = index.insert((name.get(), alt.get()));
                }
                self.index = Some(index);
            }
            None => {}
        }
    }

    pub(super) fn clear(&mut self) {
        self.names.clear();
        self.index = None;
    }
}

#[cfg(test)]
#[path = "names_tests.rs"]
mod tests;
