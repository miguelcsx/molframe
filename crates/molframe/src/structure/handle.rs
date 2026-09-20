//! The curated, owned handle onto a structure snapshot.

use molframe_core::structure::{
    AtomRef, ChainRef, ModelRef, ResidueRef, Structure as CoreStructure,
};

use super::collections::{Atoms, Chains, Models, Residues};
use super::editor::StructureEditor;
#[cfg(feature = "query")]
use super::query::QueryStructure;
#[cfg(feature = "query")]
use super::selection::Selection;
#[cfg(feature = "query")]
use molframe_core::contract::AnalysisPolicy;
#[cfg(feature = "query")]
use molframe_core::diagnostic::Findings;

/// An immutable structure snapshot.
///
/// This is a thin, `Clone`-cheap handle over the storage engine's snapshot,
/// and it is the facade's own type: what a caller can reach through it is what
/// the crate promises for 1.x, and the engine behind it is free to move. The
/// representation is private for exactly that reason — a method call, never a
/// field access, is what keeps the promise honest.
///
/// Everything a program normally wants is here: counts, the four hierarchy
/// collections, positional lookup, coordinates, [`Structure::select`] and
/// [`Structure::edit`]. When a caller needs the engine itself — a custom
/// format reader, a provider, a kernel that the facade has not yet curated —
/// [`Structure::engine`] hands over the snapshot explicitly, which is a
/// deliberate step rather than the default one.
///
/// # Examples
///
/// ```
/// use molframe::prelude::*;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # END
/// # ";
/// let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// assert_eq!(structure.chains().len(), 1);
/// assert_eq!(structure.coordinates().len(), 1);
/// assert_eq!(structure.atom_count(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Structure(CoreStructure);

impl Structure {
    /// Entry-level metadata attached to this immutable snapshot.
    #[must_use]
    pub fn metadata(&self) -> &molframe_core::structure::EntryMetadata {
        &self.0.data().entry
    }

    /// Unit cell, when the source defines one.
    #[must_use]
    pub fn cell(&self) -> Option<molframe_core::structure::UnitCell> {
        self.0.data().cell
    }

    /// Resolved bond table.
    #[must_use]
    pub fn bonds(&self) -> &molframe_core::BondTable {
        &self.0.data().bonds
    }

    /// Per-atom annotation columns.
    #[must_use]
    pub fn annotations(&self) -> &molframe_core::annotation::AtomAnnotations {
        &self.0.data().annotations
    }

    /// The positions of the first model, contiguously and without copying.
    ///
    /// This is the buffer itself, not a copy of it, which is what makes handing
    /// coordinates to a numeric array free.
    #[must_use]
    pub fn coordinates(&self) -> &[[f32; 3]] {
        self.0.positions()
    }

    /// The positions of one model, contiguously and without copying.
    ///
    /// `None` when there is no model at that position.
    #[must_use]
    pub fn model_coordinates(
        &self,
        model: molframe_core::index::ModelIndex,
    ) -> Option<&[[f32; 3]]> {
        self.0.model_positions(model)
    }

    /// A handle on every chain, in structure order.
    #[must_use]
    pub fn chains(&self) -> Chains<'_> {
        Chains::new(&self.0)
    }

    /// A handle on every residue, in structure order.
    #[must_use]
    pub fn residues(&self) -> Residues<'_> {
        Residues::new(&self.0)
    }

    /// A handle on every atom, in structure order.
    #[must_use]
    pub fn atoms(&self) -> Atoms<'_> {
        Atoms::new(&self.0)
    }

    /// A handle on every model.
    #[must_use]
    pub fn models(&self) -> Models<'_> {
        Models::new(&self.0)
    }

    /// One chain by position.
    #[must_use]
    pub fn chain_at(&self, position: usize) -> Option<ChainRef<'_>> {
        self.chains().get(position)
    }

    /// One residue by position.
    #[must_use]
    pub fn residue_at(&self, position: usize) -> Option<ResidueRef<'_>> {
        self.residues().get(position)
    }

    /// One atom by position.
    #[must_use]
    pub fn atom_at(&self, position: usize) -> Option<AtomRef<'_>> {
        self.atoms().get(position)
    }

    /// One model by position.
    #[must_use]
    pub fn model_at(&self, position: usize) -> Option<ModelRef<'_>> {
        self.models().get(position)
    }

    /// One chain by either namespace's label.
    ///
    /// The lookup is a scan of the chain table, so it is the right call for
    /// reading a structure and the wrong one for a loop over every atom.
    #[must_use]
    pub fn chain(&self, label: &str) -> Option<ChainRef<'_>> {
        self.0.data().chain_named(label)
    }

    /// The number of models.
    #[must_use]
    pub fn model_count(&self) -> usize {
        self.0.model_count()
    }

    /// The number of chains.
    #[must_use]
    pub fn chain_count(&self) -> usize {
        self.0.chain_count()
    }

    /// The number of residues.
    #[must_use]
    pub fn residue_count(&self) -> usize {
        self.0.residue_count()
    }

    /// The number of atoms in one model.
    #[must_use]
    pub fn atom_count(&self) -> u32 {
        self.0.atom_count()
    }

    /// Starts a transactional structural edit.
    #[must_use]
    pub fn edit(&self) -> StructureEditor {
        StructureEditor::new(&self.0)
    }

    /// Compiles and evaluates a selection query against this structure.
    ///
    /// # Errors
    ///
    /// Returns syntax, semantic or evaluation diagnostics.
    #[cfg(feature = "query")]
    pub fn select(&self, source: &str, policy: &AnalysisPolicy) -> Result<Selection, Findings> {
        let evaluation = QueryStructure::select_text(&self.0, source, policy)?;
        Ok(Selection::from(self.0.view_of(evaluation.selection)))
    }

    /// The storage engine's snapshot behind this handle.
    ///
    /// This is the uncurated door. What comes back is documented by
    /// `molframe::engine::core`, carries no 1.x stability promise, and is where
    /// a custom reader, a provider, or an as-yet-uncurated kernel is written
    /// against.
    #[must_use]
    pub fn engine(&self) -> &CoreStructure {
        &self.0
    }

    /// Consumes the handle and returns the engine's snapshot.
    #[must_use]
    pub fn into_engine(self) -> CoreStructure {
        self.0
    }
}

impl From<CoreStructure> for Structure {
    fn from(inner: CoreStructure) -> Self {
        Self(inner)
    }
}

impl From<Structure> for CoreStructure {
    fn from(structure: Structure) -> Self {
        structure.0
    }
}

impl std::fmt::Display for Structure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}
