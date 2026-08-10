//! Transactional copy-on-write edits over immutable structure snapshots.

use super::materialize::compact_hierarchy;
use super::{CoordinateStore, Structure, StructureData, validate};
use crate::annotation::AtomAnnotation;
use crate::bond::{BondRecord, BondTableBuilder};
use crate::chunk::ChunkBuilder;
use crate::coords::CoordinateBlock;
use crate::diagnostic::{Code, Diagnostic};
use crate::index::{AtomIndex, ChainIndex, ResidueIndex};
use crate::selection::AtomSelection;

/// A private copy-on-write snapshot that becomes visible only on commit.
///
/// Atom chunks, coordinate lanes, topology columns and the identifier arena are
/// shared with the base. Mutating one column detaches only that column.
#[derive(Debug)]
pub struct StructureEditor {
    data: StructureData,
    coordinates_changed: bool,
}

impl Structure {
    /// Starts a transactional structural edit.
    #[must_use]
    pub fn edit(&self) -> StructureEditor {
        StructureEditor {
            data: self.data().clone(),
            coordinates_changed: false,
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
        if column.len() != self.data.atom_count() {
            return Err(Diagnostic::new(Code::E3011)
                .with_context("annotation", name)
                .with_context("column rows", column.len().to_string())
                .with_context("atoms", self.data.atom_count().to_string()));
        }
        let _ = self.data.annotations.insert(name, column);
        Ok(())
    }

    /// Removes a custom annotation, returning whether it existed.
    pub fn remove_annotation(&mut self, name: &str) -> bool {
        self.data.annotations.remove(name).is_some()
    }

    /// Explicitly discards domain extensions before a topology-changing edit.
    pub fn clear_extensions(&mut self) {
        self.data.extensions.clear();
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
        self.data
            .topology
            .chains
            .rename(chain, symbol)
            .map_err(|_| Diagnostic::new(Code::E6006).with_context("chain", chain.to_string()))
    }

    /// Applies one coordinate transform to selected atoms in every dense frame.
    ///
    /// Missing coordinates remain missing. The candidate coordinate store is
    /// built privately and published to the transaction only after every
    /// transformed value has been checked.
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
        let mut coords = self.data.coords.clone();
        transform_store(&mut coords, selection, &transform)?;
        self.data.coords = coords;
        self.coordinates_changed = true;
        Ok(())
    }

    /// Deletes selected atom rows and remaps residues, coordinates and bonds.
    ///
    /// Empty residues and chains remain as declared topology. Keeping their
    /// identities avoids silently changing entity sequences or depositor
    /// numbering; their child ranges become empty and still tile atom order.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for an out-of-range atom or a ragged ensemble.
    pub fn delete_atoms(&mut self, selection: &AtomSelection) -> Result<(), Diagnostic> {
        validate_selection(selection, self.data.atom_count())?;
        if selection.is_empty() {
            return Ok(());
        }
        if matches!(self.data.coords, CoordinateStore::Ragged { .. }) {
            return Err(Diagnostic::new(Code::E6003));
        }
        require_no_extensions(&self.data)?;
        let candidate = delete_from(&self.data, selection)?;
        self.data = candidate;
        self.coordinates_changed = true;
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
        let findings = validate(&self.data);
        if findings.is_empty() {
            Ok(Structure::new(self.data))
        } else {
            Err(findings)
        }
    }
}

fn require_no_extensions(data: &StructureData) -> Result<(), Diagnostic> {
    if data.extensions.is_empty() {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E3014).with_context(
            "extensions",
            data.extensions.keys().collect::<Vec<_>>().join(","),
        ))
    }
}

fn validate_selection(selection: &AtomSelection, atom_count: u32) -> Result<(), Diagnostic> {
    if let Some(atom) = selection.iter().find(|atom| *atom >= atom_count) {
        Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()))
    } else {
        Ok(())
    }
}

fn transform_store<F>(
    store: &mut CoordinateStore,
    selection: &AtomSelection,
    transform: &F,
) -> Result<(), Diagnostic>
where
    F: Fn([f32; 3]) -> [f32; 3],
{
    match store {
        CoordinateStore::Single(block) => transform_block(block, selection, transform),
        CoordinateStore::Dense { frames } => {
            for block in frames {
                transform_block(block, selection, transform)?;
            }
            Ok(())
        }
        CoordinateStore::Ragged { .. } => Err(Diagnostic::new(Code::E6003)),
    }
}

fn transform_block<F>(
    block: &mut CoordinateBlock,
    selection: &AtomSelection,
    transform: &F,
) -> Result<(), Diagnostic>
where
    F: Fn([f32; 3]) -> [f32; 3],
{
    for atom in selection {
        let Some(position) = block.as_slice().get(atom as usize).copied() else {
            return Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()));
        };
        if !position.iter().all(|value| value.is_finite()) {
            continue;
        }
        let transformed = transform(position);
        if !transformed.iter().all(|value| value.is_finite()) {
            return Err(Diagnostic::new(Code::E6008).with_context("atom", atom.to_string()));
        }
    }
    let positions = block.as_mut_slice();
    for atom in selection {
        let Some(position) = positions.get_mut(atom as usize) else {
            return Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()));
        };
        if position.iter().all(|value| value.is_finite()) {
            *position = transform(*position);
        }
    }
    Ok(())
}

fn delete_from(data: &StructureData, deleted: &AtomSelection) -> Result<StructureData, Diagnostic> {
    let atom_count = data.atom_count();
    let mut remap = vec![u32::MAX; atom_count as usize];
    let mut next = 0u32;
    for old in 0..atom_count {
        if !deleted.contains(old) {
            remap[old as usize] = next;
            next += 1;
        }
    }

    let mut builder = ChunkBuilder::new();
    let first_positions = match data.coords.block(crate::index::ModelIndex::new(0)) {
        Some(block) => block.as_slice(),
        None => return Err(Diagnostic::new(Code::E6003)),
    };
    for chunk in data.chunks.iter() {
        builder.start_model(chunk.model());
        for old in chunk.atoms() {
            if deleted.contains(old) {
                continue;
            }
            let local = old.checked_sub(chunk.atoms().start).ok_or_else(|| {
                Diagnostic::new(Code::E3001).with_context("atom", old.to_string())
            })?;
            let position = first_positions.get(old as usize).copied();
            let Some(record) = chunk.record(local, &data.topology.residues, position) else {
                return Err(Diagnostic::new(Code::E3001).with_context("atom", old.to_string()));
            };
            builder.push(record);
        }
    }
    let (chunks, first_coords) = builder.finish();

    let mut candidate = data.clone();
    candidate.chunks = std::sync::Arc::new(chunks);
    candidate.coords = filtered_store(&data.coords, deleted, first_coords)?;
    remap_residue_ranges(&mut candidate, &remap)?;
    candidate.bonds = remap_bonds(&data.bonds, &remap);
    candidate.annotations = data
        .annotations
        .filter(|atom| !deleted.contains(atom))
        .map_err(|_| Diagnostic::new(Code::E6009))?;
    Ok(candidate)
}

fn filtered_store(
    store: &CoordinateStore,
    deleted: &AtomSelection,
    first: CoordinateBlock,
) -> Result<CoordinateStore, Diagnostic> {
    match store {
        CoordinateStore::Single(_) => Ok(CoordinateStore::Single(first)),
        CoordinateStore::Dense { frames } => {
            let mut filtered = Vec::with_capacity(frames.len());
            for (index, frame) in frames.iter().enumerate() {
                if index == 0 {
                    filtered.push(first.clone());
                } else {
                    filtered.push(
                        frame
                            .as_slice()
                            .iter()
                            .copied()
                            .enumerate()
                            .filter_map(|(atom, position)| {
                                let atom = u32::try_from(atom).ok()?;
                                (!deleted.contains(atom)).then_some(position)
                            })
                            .collect(),
                    );
                }
            }
            Ok(CoordinateStore::Dense { frames: filtered })
        }
        CoordinateStore::Ragged { .. } => Err(Diagnostic::new(Code::E6003)),
    }
}

fn remap_residue_ranges(data: &mut StructureData, remap: &[u32]) -> Result<(), Diagnostic> {
    let mut first = 0u32;
    for position in
        0..u32::try_from(data.topology.residues.len()).map_err(|_| Diagnostic::new(Code::E6009))?
    {
        let residue = ResidueIndex::new(position);
        let Some(old) = data.topology.residues.atoms(residue) else {
            continue;
        };
        let kept = old
            .filter(|atom| {
                remap
                    .get(*atom as usize)
                    .is_some_and(|new| *new != u32::MAX)
            })
            .count();
        let kept = u32::try_from(kept).map_err(|_| Diagnostic::new(Code::E6009))?;
        data.topology
            .residues
            .set_atoms(
                residue,
                first
                    ..first
                        .checked_add(kept)
                        .ok_or_else(|| Diagnostic::new(Code::E6009))?,
            )
            .map_err(|_| Diagnostic::new(Code::E3001))?;
        first = first
            .checked_add(kept)
            .ok_or_else(|| Diagnostic::new(Code::E6009))?;
    }
    Ok(())
}

fn remap_bonds(table: &crate::bond::BondTable, remap: &[u32]) -> crate::bond::BondTable {
    if !table.is_available() {
        return crate::bond::BondTable::default();
    }
    let mut builder = BondTableBuilder::new();
    for bond in table.iter() {
        let Some(atom_a) = remapped(bond.atom_a, remap) else {
            continue;
        };
        let Some(atom_b) = remapped(bond.atom_b, remap) else {
            continue;
        };
        builder.push(BondRecord {
            atom_a,
            atom_b,
            order: bond.order,
            provenance: bond.provenance,
        });
    }
    builder.finish()
}

fn remapped(atom: AtomIndex, remap: &[u32]) -> Option<AtomIndex> {
    let mapped = *remap.get(atom.as_usize())?;
    (mapped != u32::MAX).then(|| AtomIndex::new(mapped))
}

#[cfg(test)]
#[path = "overlay_tests.rs"]
mod tests;
