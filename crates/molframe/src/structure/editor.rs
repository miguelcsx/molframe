//! The curated entry point for transactional structural edits.

use molframe_core::ExecutionContext;
use molframe_core::MemoryBudgetError;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::index::ChainIndex;
use molframe_core::structure::{
    CoordinateEditor, Structure as CoreStructure, StructureEditor as CoreEditor,
};

use super::handle::Structure;
use super::selection::Selection;

/// A transactional structural edit, scoped to one snapshot.
///
/// Topology edits (deleting atoms, renaming a chain) are staged here and
/// published together by [`StructureEditor::finish`]. A scoped, memory
/// -budgeted coordinate rewrite is a separate transaction reached through
/// [`StructureEditor::coordinates`] — coordinates and topology are edited
/// through different engine mechanisms today, so the facade exposes both
/// rather than pretending they are one.
///
/// # Examples
///
/// ```
/// use molframe::prelude::*;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
/// # END
/// # ";
/// let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let mut editor = structure.edit();
/// editor
///     .rename_chain(molframe::ChainIndex::new(0), "B")
///     .map_err(Findings::from)?;
/// let renamed = editor.finish().map_err(Findings::from)?;
/// assert_eq!(renamed.chain_at(0).and_then(|chain| chain.label()), Some("B"));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct StructureEditor {
    base: CoreStructure,
    overlay: CoreEditor,
}

impl StructureEditor {
    pub(super) fn new(base: &CoreStructure) -> Self {
        Self {
            base: base.clone(),
            overlay: base.edit(),
        }
    }

    /// Stages the deletion of every atom `selection` covers.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an out-of-range atom or a ragged ensemble.
    pub fn delete(&mut self, selection: &Selection) -> Result<(), Diagnostic> {
        self.overlay.delete_atoms(selection.atom_selection())
    }

    /// Stages a rename of one chain in both identifier namespaces.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the chain does not exist, the new label is
    /// empty, or the bounded identifier dictionary cannot accept it.
    pub fn rename_chain(&mut self, chain: ChainIndex, label: &str) -> Result<(), Diagnostic> {
        self.overlay.rename_chain(chain, label)
    }

    /// Starts a scoped, memory-budgeted coordinate rewrite of the original
    /// snapshot this editor was opened from.
    ///
    /// # Errors
    ///
    /// Returns before publishing an editor when all coordinate copies cannot
    /// fit.
    pub fn coordinates(
        &self,
        context: &ExecutionContext,
    ) -> Result<CoordinateEditor, MemoryBudgetError> {
        self.base.edit_coordinates(context)
    }

    /// Validates the staged topology edits and publishes a new snapshot.
    ///
    /// # Errors
    ///
    /// Returns every violated structural invariant. The original structure is
    /// unaffected whether this succeeds or fails.
    pub fn finish(self) -> Result<Structure, Vec<Diagnostic>> {
        self.overlay.commit().map(Structure::from)
    }
}
