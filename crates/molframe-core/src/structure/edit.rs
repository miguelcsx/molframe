//! Transactional coordinate editing over immutable snapshots.

use super::{CoordinateStore, Structure, StructureData, validate};
use crate::diagnostic::{Code, Diagnostic};
use crate::index::ModelIndex;
use crate::{ExecutionContext, MemoryBudgetError};

/// A private coordinate copy that becomes a new structure only on commit.
///
/// The original structure remains immutable and every view into it stays valid.
/// Dropping an editor discards the transaction.
#[derive(Debug)]
pub struct CoordinateEditor {
    base: Structure,
    coords: CoordinateStore,
}

impl Structure {
    /// Starts a scoped coordinate edit under the shared execution account.
    ///
    /// Coordinate storage is copied for the new snapshot; topology and atom
    /// annotations remain untouched. Call [`CoordinateEditor::commit`] to
    /// validate and publish the new snapshot.
    /// Shared topology remains borrowed. Direct coordinate frames detach under
    /// the supplied account, and their reservations follow the copied backing
    /// buffers through committed snapshots and clones.
    ///
    /// # Errors
    ///
    /// Returns before publishing an editor when all coordinate copies cannot fit.
    pub fn edit_coordinates(
        &self,
        context: &ExecutionContext,
    ) -> Result<CoordinateEditor, MemoryBudgetError> {
        let mut coords = self.data().coords.clone();
        coords.make_unique_in(context)?;
        Ok(CoordinateEditor {
            base: self.clone(),
            coords,
        })
    }
}

impl CoordinateEditor {
    /// Mutable positions of one dense model.
    ///
    /// Ragged models are independent structures and must be edited through the
    /// corresponding value returned by [`Structure::ragged_models`].
    pub fn positions_mut(&mut self, model: ModelIndex) -> Option<&mut [[f32; 3]]> {
        match &mut self.coords {
            CoordinateStore::Single(block) if model.get() == 0 => Some(block.as_mut_slice()),
            CoordinateStore::Dense { frames } => frames
                .get_mut(model.as_usize())
                .map(crate::coords::CoordinateBlock::as_mut_slice),
            CoordinateStore::Single(_) | CoordinateStore::Ragged { .. } => None,
        }
    }

    /// Mutable positions of one dense model, with a diagnostic on failure.
    ///
    /// This is the fallible form for callers that must preserve molframe's
    /// diagnostic contract across an FFI boundary.
    ///
    /// # Errors
    ///
    /// Returns `E6003` when `model` is not a dense frame of this structure.
    pub fn try_positions_mut(&mut self, model: ModelIndex) -> Result<&mut [[f32; 3]], Diagnostic> {
        self.positions_mut(model)
            .ok_or_else(|| Diagnostic::new(Code::E6003))
    }

    /// Validates the transaction and returns a new immutable snapshot.
    ///
    /// # Errors
    ///
    /// Returns every violated structural invariant. The original snapshot is
    /// unaffected whether the commit succeeds or fails.
    pub fn commit(self) -> Result<Structure, Vec<crate::diagnostic::Diagnostic>> {
        let mut data: StructureData = self.base.data().clone();
        data.coords = self.coords;
        finish(data)
    }

    /// Builds a validated snapshot while retaining the editor's private buffer.
    ///
    /// This is useful at an FFI boundary where a foreign array may outlive its
    /// editing scope. The published structure receives its own coordinate
    /// buffer, so a retained foreign view cannot mutate it afterwards.
    ///
    /// # Errors
    ///
    /// Returns every violated structural invariant.
    pub fn snapshot(&self) -> Result<Structure, Vec<crate::diagnostic::Diagnostic>> {
        let mut data: StructureData = self.base.data().clone();
        data.coords = self.coords.clone();
        finish(data)
    }
}

fn finish(mut data: StructureData) -> Result<Structure, Vec<crate::diagnostic::Diagnostic>> {
    let Some(generation) = data.generation.next() else {
        return Err(vec![
            Diagnostic::new(Code::E6003).with_context("coordinate_generation", "exhausted"),
        ]);
    };
    data.generation = generation;
    let findings = validate(&data);
    if findings.is_empty() {
        Ok(Structure::new(data))
    } else {
        Err(findings)
    }
}

#[cfg(test)]
#[path = "edit_tests.rs"]
mod tests;
