//! Hierarchy expansion for `byres` and `same`.

use pdbiox_core::index::{ChainIndex, EntityIndex, ResidueIndex};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use std::collections::BTreeSet;

pub(super) fn same_residue(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
) -> AtomSelection {
    let residues: BTreeSet<ResidueIndex> = selected
        .iter()
        .filter_map(|atom| structure.atom(pdbiox_core::index::AtomIndex::new(atom)))
        .filter_map(pdbiox_core::structure::AtomRef::residue)
        .map(pdbiox_core::structure::ResidueRef::index)
        .collect();
    select_atoms(structure, universe, |residue, _, _| {
        residues.contains(&residue)
    })
}

pub(super) fn same_chain(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
) -> AtomSelection {
    let selected_residues: BTreeSet<ResidueIndex> = selected
        .iter()
        .filter_map(|atom| structure.atom(pdbiox_core::index::AtomIndex::new(atom)))
        .filter_map(pdbiox_core::structure::AtomRef::residue)
        .map(pdbiox_core::structure::ResidueRef::index)
        .collect();
    let chains: BTreeSet<ChainIndex> = structure
        .data()
        .chains()
        .filter(|chain| {
            chain
                .residues()
                .any(|residue| selected_residues.contains(&residue.index()))
        })
        .map(pdbiox_core::structure::ChainRef::index)
        .collect();
    select_atoms(structure, universe, |_, chain, _| chains.contains(&chain))
}

pub(super) fn same_entity(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
) -> AtomSelection {
    let selected_residues: BTreeSet<ResidueIndex> = selected
        .iter()
        .filter_map(|atom| structure.atom(pdbiox_core::index::AtomIndex::new(atom)))
        .filter_map(pdbiox_core::structure::AtomRef::residue)
        .map(pdbiox_core::structure::ResidueRef::index)
        .collect();
    let entities: BTreeSet<EntityIndex> = structure
        .data()
        .chains()
        .filter(|chain| {
            chain
                .residues()
                .any(|residue| selected_residues.contains(&residue.index()))
        })
        .filter_map(pdbiox_core::structure::ChainRef::entity)
        .collect();
    select_atoms(structure, universe, |_, _, entity| {
        entities.contains(&entity)
    })
}

fn select_atoms(
    structure: &Structure,
    universe: &AtomSelection,
    mut accepts: impl FnMut(ResidueIndex, ChainIndex, EntityIndex) -> bool,
) -> AtomSelection {
    let mut positions = Vec::new();
    for chain in structure.data().chains() {
        let Some(entity) = chain.entity() else {
            continue;
        };
        for residue in chain.residues() {
            if !accepts(residue.index(), chain.index(), entity) {
                continue;
            }
            positions.extend(
                residue
                    .atoms()
                    .map(pdbiox_core::structure::AtomRef::index)
                    .map(pdbiox_core::index::AtomIndex::get)
                    .filter(|atom| universe.contains(*atom)),
            );
        }
    }
    AtomSelection::from_sorted(positions)
}
