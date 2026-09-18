//! Internal coordinate, hierarchy, and bond remapping helpers for edits.

use super::super::{CoordinateStore, StructureData};
use super::DELETED_ATOM;
use crate::bond::{BondRecord, BondTableBuilder};
use crate::chunk::ChunkBuilder;
use crate::coords::CoordinateBlock;
use crate::diagnostic::{Code, Diagnostic};
use crate::index::{AtomIndex, ResidueIndex};
use crate::selection::AtomSelection;

pub(super) fn require_no_extensions(data: &StructureData) -> Result<(), Diagnostic> {
    if data.extensions.is_empty() {
        return Ok(());
    }
    let mut names = String::new();
    for name in data.extensions.keys() {
        if !names.is_empty() {
            names.push(',');
        }
        names.push_str(name);
    }
    Err(Diagnostic::new(Code::E3014).with_context("extensions", names))
}

pub(super) fn validate_selection(
    selection: &AtomSelection,
    atom_count: u32,
) -> Result<(), Diagnostic> {
    if let Some(atom) = selection.iter().find(|atom| *atom >= atom_count) {
        Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()))
    } else {
        Ok(())
    }
}

pub(super) fn transform_store<F>(
    store: &mut CoordinateStore,
    selection: &AtomSelection,
    transform: &F,
) -> Result<bool, Diagnostic>
where
    F: Fn([f32; 3]) -> [f32; 3],
{
    match store {
        CoordinateStore::Single(block) => transform_block(block, selection, transform),
        CoordinateStore::Dense { frames } => {
            let mut changed = false;
            for block in frames {
                if transform_block(block, selection, transform)? {
                    changed = true;
                }
            }
            Ok(changed)
        }
        CoordinateStore::Ragged { .. } => Err(Diagnostic::new(Code::E6003)),
    }
}

fn transform_block<F>(
    block: &mut CoordinateBlock,
    selection: &AtomSelection,
    transform: &F,
) -> Result<bool, Diagnostic>
where
    F: Fn([f32; 3]) -> [f32; 3],
{
    let positions = block.as_mut_slice();
    let mut changed = false;
    for atom in selection {
        let Some(position) = positions.get_mut(atom as usize) else {
            return Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()));
        };
        if !position.iter().all(|value| value.is_finite()) {
            continue;
        }
        let transformed = transform(*position);
        if !transformed.iter().all(|value| value.is_finite()) {
            return Err(Diagnostic::new(Code::E6008).with_context("atom", atom.to_string()));
        }
        *position = transformed;
        changed = true;
    }
    Ok(changed)
}

fn build_atom_remap(atom_count: u32, deleted: &AtomSelection) -> Result<Vec<u32>, Diagnostic> {
    let len = usize::try_from(atom_count).map_err(|_| Diagnostic::new(Code::E6009))?;
    let mut remap = vec![0; len];
    for atom in deleted {
        let Some(slot) = remap.get_mut(atom as usize) else {
            return Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()));
        };
        *slot = DELETED_ATOM;
    }
    let mut next = 0u32;
    for mapped in &mut remap {
        if *mapped == DELETED_ATOM {
            continue;
        }
        *mapped = next;
        next = next
            .checked_add(1)
            .ok_or_else(|| Diagnostic::new(Code::E6009))?;
    }
    Ok(remap)
}

pub(super) fn delete_from(
    data: &StructureData,
    deleted: &AtomSelection,
) -> Result<StructureData, Diagnostic> {
    let remap = build_atom_remap(data.atom_count(), deleted)?;
    let first_positions = match data.coords.block(crate::index::ModelIndex::new(0)) {
        Some(block) => block.as_slice(),
        None => return Err(Diagnostic::new(Code::E6003)),
    };
    let mut builder = ChunkBuilder::new();
    for chunk in data.chunks.iter() {
        builder.start_model(chunk.model());
        for old in chunk.atoms() {
            let mapped = remap.get(old as usize).copied().ok_or_else(|| {
                Diagnostic::new(Code::E3001).with_context("atom", old.to_string())
            })?;
            if mapped == DELETED_ATOM {
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
    candidate.coords = filtered_store(&data.coords, &remap, first_coords)?;
    remap_residue_ranges(&mut candidate, &remap)?;
    candidate.bonds = remap_bonds(&data.bonds, &remap);
    candidate.annotations = data
        .annotations
        .filter(|atom| {
            remap
                .get(atom as usize)
                .is_some_and(|mapped| *mapped != DELETED_ATOM)
        })
        .map_err(|_| Diagnostic::new(Code::E6009))?;
    Ok(candidate)
}

fn filtered_store(
    store: &CoordinateStore,
    remap: &[u32],
    first: CoordinateBlock,
) -> Result<CoordinateStore, Diagnostic> {
    match store {
        CoordinateStore::Single(block) => {
            validate_coordinate_rows(block, remap.len())?;
            Ok(CoordinateStore::Single(first))
        }
        CoordinateStore::Dense { frames } => {
            if frames.is_empty() {
                return Err(Diagnostic::new(Code::E6003));
            }
            for frame in frames {
                validate_coordinate_rows(frame, remap.len())?;
            }
            let mut filtered = Vec::with_capacity(frames.len());
            filtered.push(first);
            for frame in frames.iter().skip(1) {
                filtered.push(filter_coordinate_block(frame, remap));
            }
            Ok(CoordinateStore::Dense { frames: filtered })
        }
        CoordinateStore::Ragged { .. } => Err(Diagnostic::new(Code::E6003)),
    }
}

fn validate_coordinate_rows(block: &CoordinateBlock, expected: usize) -> Result<(), Diagnostic> {
    if block.as_slice().len() == expected {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E6003)
            .with_context("coordinate rows", block.as_slice().len().to_string())
            .with_context("atoms", expected.to_string()))
    }
}

fn filter_coordinate_block(frame: &CoordinateBlock, remap: &[u32]) -> CoordinateBlock {
    frame
        .as_slice()
        .iter()
        .copied()
        .zip(remap.iter().copied())
        .filter_map(|(position, mapped)| (mapped != DELETED_ATOM).then_some(position))
        .collect()
}

fn remap_residue_ranges(data: &mut StructureData, remap: &[u32]) -> Result<(), Diagnostic> {
    let residue_count =
        u32::try_from(data.topology.residues.len()).map_err(|_| Diagnostic::new(Code::E6009))?;
    let mut first = 0u32;
    for position in 0..residue_count {
        let residue = ResidueIndex::new(position);
        let Some(old) = data.topology.residues.atoms(residue) else {
            continue;
        };
        let start = usize::try_from(old.start).map_err(|_| Diagnostic::new(Code::E3001))?;
        let end = usize::try_from(old.end).map_err(|_| Diagnostic::new(Code::E3001))?;
        let old_mapping = remap
            .get(start..end)
            .ok_or_else(|| Diagnostic::new(Code::E3001))?;
        let kept = old_mapping
            .iter()
            .filter(|mapped| **mapped != DELETED_ATOM)
            .count();
        let kept = u32::try_from(kept).map_err(|_| Diagnostic::new(Code::E6009))?;
        let end = first
            .checked_add(kept)
            .ok_or_else(|| Diagnostic::new(Code::E6009))?;
        data.topology
            .residues
            .set_atoms(residue, first..end)
            .map_err(|_| Diagnostic::new(Code::E3001))?;
        first = end;
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
    match remap.get(atom.as_usize()).copied() {
        Some(DELETED_ATOM) | None => None,
        Some(mapped) => Some(AtomIndex::new(mapped)),
    }
}
