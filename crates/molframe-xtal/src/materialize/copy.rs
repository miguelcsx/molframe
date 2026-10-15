//! Counting, copying and remapping for assembly materialisation.
//!
//! Complexity: O(instances) for the counts; O(atoms) per annotation or coordinate
//! column copied; O(1) per atom or bond remap.

use super::errors::{capacity, invariant};
use super::{CopyGroup, CopySpan, TransformKey};
use crate::AssemblyView;
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::bond::BondRecord;
use molframe_core::chunk::AtomRecord;
use molframe_core::optional::OptionalSymbol;
use molframe_core::symbol::SymbolId;
use molframe_core::{AtomIndex, ChainIndex, Code, Diagnostic, Structure};
use molframe_geom::Rigid;
use std::collections::BTreeMap;

/// Computes the expected expanded atom count with overflow checking.
pub(super) fn expanded_atom_count(view: &AssemblyView) -> Result<usize, Diagnostic> {
    view.instances.iter().try_fold(0usize, |total, instance| {
        total
            .checked_add(instance.atoms.len())
            .ok_or_else(|| capacity("atoms"))
    })
}

/// Reconstructs one source atom record using a monotonic chunk cursor.
///
/// Assembly instance atoms are copied in source order, allowing chunk lookup to
/// advance linearly instead of binary-searching the chunk table for every atom.
pub(super) fn atom_record(
    source: &Structure,
    atom: AtomIndex,
    chunk_cursor: &mut usize,
) -> Result<AtomRecord, Diagnostic> {
    let chunks = &source.data().chunks;

    loop {
        let chunk = chunks.get(*chunk_cursor).ok_or_else(invariant)?;
        let atoms = chunk.atoms();

        if atom.get() < atoms.start {
            return Err(invariant());
        }

        if atom.get() < atoms.end {
            let local = atom.get().checked_sub(atoms.start).ok_or_else(invariant)?;
            let position = source.positions().get(atom.as_usize()).copied();

            return chunk
                .record(local, &source.data().topology.residues, position)
                .ok_or_else(invariant);
        }

        *chunk_cursor = chunk_cursor
            .checked_add(1)
            .ok_or_else(|| capacity("chunks"))?;
    }
}

/// Copies one annotation column according to the output-to-source atom mapping.
pub(super) fn copy_annotation(
    annotation: &AtomAnnotation,
    source_atoms: &[u32],
) -> Result<AtomAnnotation, Diagnostic> {
    match annotation {
        AtomAnnotation::Boolean(column) => {
            copy_column(column, source_atoms).map(AtomAnnotation::Boolean)
        }
        AtomAnnotation::Integer(column) => {
            copy_column(column, source_atoms).map(AtomAnnotation::Integer)
        }
        AtomAnnotation::Real(column) => copy_column(column, source_atoms).map(AtomAnnotation::Real),
        AtomAnnotation::Symbol(column) => {
            copy_column(column, source_atoms).map(AtomAnnotation::Symbol)
        }
        _ => {
            Err(Diagnostic::new(Code::E3013)
                .with_context("annotation", "unsupported physical type"))
        }
    }
}

/// Copies entries from one physical annotation column in output atom order.
pub(super) fn copy_column<T: Copy>(
    column: &AnnotationColumn<T>,
    source_atoms: &[u32],
) -> Result<AnnotationColumn<T>, Diagnostic> {
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(source_atoms.len())
        .map_err(|_| capacity("annotation rows"))?;

    for &atom in source_atoms {
        entries.push(column.get(atom).ok_or_else(invariant)?);
    }

    AnnotationColumn::from_entries(entries)
        .map_err(|error| invariant().with_context("cause", error.to_string()))
}

/// Groups assembly copies sharing an exactly identical rigid transform.
pub(super) fn copy_groups(copies: &[CopySpan]) -> Vec<CopyGroup> {
    let mut grouped = BTreeMap::<TransformKey, CopyGroup>::new();

    for (position, copy) in copies.iter().enumerate() {
        grouped
            .entry(transform_key(copy.transform))
            .or_default()
            .entry(copy.source_chain.get())
            .or_default()
            .push(position);
    }

    grouped.into_values().collect()
}

/// Finds transform groups containing both source chains of a bond.
pub(super) fn groups_for_pair(groups: &[CopyGroup], chain_a: u32, chain_b: u32) -> Vec<usize> {
    groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            (group.contains_key(&chain_a) && group.contains_key(&chain_b)).then_some(index)
        })
        .collect()
}

/// Produces an orientation-independent cache key for a source-chain pair.
pub(super) fn chain_pair_key(left: ChainIndex, right: ChainIndex) -> (u32, u32) {
    let left = left.get();
    let right = right.get();

    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

/// Produces an exact deterministic key for a rigid transform.
pub(super) fn transform_key(transform: Rigid) -> TransformKey {
    let mut key = [0_u64; 12];

    for row in 0..3 {
        for column in 0..3 {
            key[row * 3 + column] = transform.rotation[row][column].to_bits();
        }

        key[9 + row] = transform.translation[row].to_bits();
    }

    key
}

/// Resolves the source chain owning one atom.
pub(super) fn atom_chain(source: &Structure, atom: AtomIndex) -> Result<ChainIndex, Diagnostic> {
    let residue = source
        .data()
        .topology
        .residues
        .containing(atom.get())
        .ok_or_else(invariant)?;

    source
        .data()
        .topology
        .chains
        .containing(residue.get())
        .ok_or_else(invariant)
}

/// Remaps both endpoints of one source bond into materialized copies.
///
/// The assembly shares the source structure's dictionary, so the connection's
/// own classification stays valid and is copied with the rest of the record.
pub(super) fn remap_bond(
    bond: BondRecord,
    left: &CopySpan,
    right: &CopySpan,
) -> Result<BondRecord, Diagnostic> {
    Ok(BondRecord {
        atom_a: remap_atom(bond.atom_a, left)?,
        atom_b: remap_atom(bond.atom_b, right)?,
        ..bond
    })
}

/// Remaps one source atom into the output span of a copied chain.
pub(super) fn remap_atom(atom: AtomIndex, copy: &CopySpan) -> Result<AtomIndex, Diagnostic> {
    if !copy.source_atoms.contains(&atom.get()) {
        return Err(invariant());
    }

    let local = atom
        .get()
        .checked_sub(copy.source_atoms.start)
        .ok_or_else(invariant)?;

    copy.output_start
        .checked_add(local)
        .map(AtomIndex::new)
        .ok_or_else(|| capacity("atoms"))
}

/// Converts an optional symbol into its compact topology representation.
pub(super) fn optional(symbol: Option<SymbolId>) -> OptionalSymbol {
    match symbol {
        Some(symbol) => OptionalSymbol::some(symbol),
        None => OptionalSymbol::NONE,
    }
}

pub(super) fn expanded_residue_count(view: &AssemblyView) -> Result<usize, Diagnostic> {
    view.instances.iter().try_fold(0usize, |total, instance| {
        let residues = view
            .source()
            .data()
            .topology
            .chains
            .residues(instance.source_chain)
            .ok_or_else(invariant)?;
        total
            .checked_add(residues.len())
            .ok_or_else(|| capacity("residues"))
    })
}

/// Allocates an empty vector with fallible exact capacity for the named resource.
pub(super) fn reserved<T>(count: usize, name: &str) -> Result<Vec<T>, Diagnostic> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| capacity(name))?;
    Ok(values)
}

pub(super) fn secondary_storage<T>(count: usize, present: bool) -> Result<Vec<T>, Diagnostic> {
    if present {
        reserved(count, "secondary rows")
    } else {
        Ok(Vec::new())
    }
}
