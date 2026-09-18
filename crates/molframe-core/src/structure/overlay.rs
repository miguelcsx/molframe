//! Transactional copy-on-write edits over immutable structure snapshots.

mod helpers;

use super::materialize::compact_hierarchy;
use super::{Structure, StructureData, validate, validate_coordinate_edit};
use crate::annotation::AtomAnnotation;
#[cfg(test)]
use crate::bond::{BondRecord, BondTableBuilder};
use crate::diagnostic::{Code, Diagnostic};
use crate::index::ChainIndex;
#[cfg(test)]
use crate::index::{AtomIndex, ResidueIndex};
use crate::selection::AtomSelection;
use helpers::{delete_from, require_no_extensions, transform_store, validate_selection};

const DELETED_ATOM: u32 = u32::MAX;

/// A private copy-on-write snapshot that becomes visible only on commit.
#[derive(Debug)]
pub struct StructureEditor {
    data: StructureData,
    coordinates_changed: bool,
    requires_full_validation: bool,
}

impl Structure {
    /// Starts a transactional structural edit.
    #[must_use]
    pub fn edit(&self) -> StructureEditor {
        StructureEditor {
            data: self.data().clone(),
            coordinates_changed: false,
            requires_full_validation: false,
        }
    }

    /// Copies a selected atom set into a compact, independent structure.
    ///
    /// Empty residues, chains and unreferenced entities are removed. Domain
    /// extensions are deliberately discarded because their topology-aligned
    /// contents cannot remain valid after materialisation.
    ///
    /// # Errors
    ///
    /// Returns diagnostics for an out-of-range atom, a ragged ensemble, or an
    /// invalid source hierarchy.
    pub fn materialize(&self, selection: &AtomSelection) -> Result<Structure, Vec<Diagnostic>> {
        if let Err(error) = validate_selection(selection, self.data().atom_count()) {
            return Err(vec![error]);
        }

        let deleted = AtomSelection::All(self.atom_count()).difference(selection);
        let mut editor = self.edit();
        editor.clear_extensions();
        if let Err(error) = editor.delete_atoms(&deleted) {
            return Err(vec![error]);
        }
        if let Err(error) = compact_hierarchy(&mut editor.data) {
            return Err(vec![error]);
        }
        editor.commit()
    }
}

impl StructureEditor {
    /// Adds or replaces a typed per-atom annotation column.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the name is empty or the column length differs
    /// from the structure's atom count.
    pub fn set_annotation(&mut self, name: &str, column: AtomAnnotation) -> Result<(), Diagnostic> {
        if name.is_empty() {
            return Err(Diagnostic::new(Code::E6007));
        }
        let atom_count = self.data.atom_count();
        if column.len() != atom_count {
            return Err(Diagnostic::new(Code::E3011)
                .with_context("annotation", name)
                .with_context("column rows", column.len().to_string())
                .with_context("atoms", atom_count.to_string()));
        }
        let _ = self.data.annotations.insert(name, column);
        self.requires_full_validation = true;
        Ok(())
    }

    /// Removes a custom annotation, returning whether it existed.
    pub fn remove_annotation(&mut self, name: &str) -> bool {
        let removed = self.data.annotations.remove(name).is_some();
        self.requires_full_validation |= removed;
        removed
    }

    /// Explicitly discards domain extensions before a topology-changing edit.
    pub fn clear_extensions(&mut self) {
        self.data.extensions.clear();
        self.requires_full_validation = true;
    }

    /// Renames a chain in both identifier namespaces.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the chain does not exist, the new label is
    /// empty, or the bounded identifier dictionary cannot accept it.
    pub fn rename_chain(&mut self, chain: ChainIndex, label: &str) -> Result<(), Diagnostic> {
        if self.data.topology.chains.label_asym_id(chain).is_none() {
            return Err(Diagnostic::new(Code::E6006).with_context("chain", chain.to_string()));
        }
        if label.is_empty() {
            return Err(Diagnostic::new(Code::E6007));
        }
        require_no_extensions(&self.data)?;
        let symbol = self.data.dictionary.intern(label).map_err(|_| {
            Diagnostic::new(Code::E1901)
                .with_context("limit", "dictionary entries")
                .with_context("identifier", label)
        })?;
        let result = self
            .data
            .topology
            .chains
            .rename(chain, symbol)
            .map_err(|_| Diagnostic::new(Code::E6006).with_context("chain", chain.to_string()));
        if result.is_ok() {
            self.requires_full_validation = true;
        }
        result
    }

    /// Applies one coordinate transform to selected atoms in every dense frame.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an out-of-range atom, a non-finite result, or a
    /// ragged ensemble whose models must be edited as independent snapshots.
    pub fn transform<F>(
        &mut self,
        selection: &AtomSelection,
        transform: F,
    ) -> Result<(), Diagnostic>
    where
        F: Fn([f32; 3]) -> [f32; 3],
    {
        validate_selection(selection, self.data.atom_count())?;
        if selection.is_empty() {
            return Ok(());
        }
        let mut coords = self.data.coords.clone();
        let changed = transform_store(&mut coords, selection, &transform)?;
        if changed {
            self.data.coords = coords;
            self.coordinates_changed = true;
        }
        Ok(())
    }

    /// Deletes selected atom rows and remaps residues, coordinates and bonds.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an out-of-range atom or a ragged ensemble.
    pub fn delete_atoms(&mut self, selection: &AtomSelection) -> Result<(), Diagnostic> {
        validate_selection(selection, self.data.atom_count())?;
        if selection.is_empty() {
            return Ok(());
        }
        if matches!(&self.data.coords, super::CoordinateStore::Ragged { .. }) {
            return Err(Diagnostic::new(Code::E6003));
        }
        require_no_extensions(&self.data)?;
        self.data = delete_from(&self.data, selection)?;
        self.coordinates_changed = true;
        self.requires_full_validation = true;
        Ok(())
    }

    /// Validates and publishes the edited snapshot.
    ///
    /// # Errors
    ///
    /// Returns every violated structural invariant. The original structure is
    /// unaffected whether commit succeeds or fails.
    pub fn commit(mut self) -> Result<Structure, Vec<Diagnostic>> {
        if self.coordinates_changed {
            let Some(generation) = self.data.generation.next() else {
                return Err(vec![
                    Diagnostic::new(Code::E6003).with_context("coordinate_generation", "exhausted"),
                ]);
            };
            self.data.generation = generation;
        }
        let findings = if self.requires_full_validation {
            validate(&self.data)
        } else {
            validate_coordinate_edit(&self.data)
        };
        if findings.is_empty() {
            Ok(Structure::new(self.data))
        } else {
            Err(findings)
        }
    }
}

#[cfg(test)]
#[path = "overlay_tests.rs"]
mod tests;
