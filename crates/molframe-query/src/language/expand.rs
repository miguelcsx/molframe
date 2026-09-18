//! Hierarchy expansion for `byres` and `same`.

use molframe_core::index::{AtomIndex, ChainIndex, EntityIndex, ResidueIndex};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{ResidueRef, Structure};
use std::collections::HashSet;

/// Expands `selected` to complete residues inside `universe`.
///
/// Selected residue discovery is `O(S)` expected time for `S` selected atoms.
/// The hierarchy expansion is `O(R + M)` for visited residues and atoms, with
/// `O(K)` auxiliary space for `K` distinct selected residues.
pub(crate) fn same_residue(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
) -> AtomSelection {
    if selected.is_empty() {
        return AtomSelection::Empty;
    }

    let residues = selected_residues(structure, selected);

    if residues.is_empty() {
        return AtomSelection::Empty;
    }

    select_residues(structure, universe, &residues)
}

/// Expands `selected` to complete chains inside `universe`.
///
/// Chain identifiers are derived directly from selected atoms instead of
/// rescanning every residue to discover their containing chains. Discovery is
/// `O(S)` expected time and expansion performs one hierarchy traversal.
pub(crate) fn same_chain(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
) -> AtomSelection {
    if selected.is_empty() {
        return AtomSelection::Empty;
    }

    let chains = selected_chains(structure, selected);

    if chains.is_empty() {
        return AtomSelection::Empty;
    }

    select_chains(structure, universe, &chains)
}

/// Expands `selected` to all atoms belonging to the same entities.
///
/// Entity identifiers are derived directly from the selected atoms through
/// residue-to-chain topology lookup. Discovery is `O(S)` expected time and
/// expansion performs one hierarchy traversal.
pub(crate) fn same_entity(
    structure: &Structure,
    universe: &AtomSelection,
    selected: &AtomSelection,
) -> AtomSelection {
    if selected.is_empty() {
        return AtomSelection::Empty;
    }

    let entities = selected_entities(structure, selected);

    if entities.is_empty() {
        return AtomSelection::Empty;
    }

    select_entities(structure, universe, &entities)
}

/// Collects distinct residue indices referenced by `selected`.
///
/// Invalid or hierarchy-less atom indices are ignored, preserving the original
/// behavior. Expected runtime is `O(S)` with `O(K)` space.
fn selected_residues(structure: &Structure, selected: &AtomSelection) -> HashSet<ResidueIndex> {
    let mut residues = HashSet::new();

    for atom in selected {
        let Some(residue) = structure
            .atom(AtomIndex::new(atom))
            .and_then(molframe_core::structure::AtomRef::residue)
        else {
            continue;
        };

        residues.insert(residue.index());
    }

    residues
}

/// Collects distinct containing-chain indices referenced by `selected`.
///
/// The topology's direct residue-to-chain mapping removes the previous
/// full-chain/residue discovery scan. Expected runtime is `O(S)`.
fn selected_chains(structure: &Structure, selected: &AtomSelection) -> HashSet<ChainIndex> {
    let mut chains = HashSet::new();

    for atom in selected {
        let Some(residue) = structure
            .atom(AtomIndex::new(atom))
            .and_then(molframe_core::structure::AtomRef::residue)
        else {
            continue;
        };

        let Some(chain) = structure
            .data()
            .topology
            .chains
            .containing(residue.index().get())
        else {
            continue;
        };

        chains.insert(chain);
    }

    chains
}

/// Collects distinct entity indices referenced by `selected`.
///
/// Discovery follows atom → residue → chain → entity directly and therefore
/// runs in expected `O(S)` time with `O(K)` auxiliary space.
fn selected_entities(structure: &Structure, selected: &AtomSelection) -> HashSet<EntityIndex> {
    let mut entities = HashSet::new();

    for atom in selected {
        let Some(residue) = structure
            .atom(AtomIndex::new(atom))
            .and_then(molframe_core::structure::AtomRef::residue)
        else {
            continue;
        };

        let Some(chain_index) = structure
            .data()
            .topology
            .chains
            .containing(residue.index().get())
        else {
            continue;
        };

        let Some(chain) = structure.chain(chain_index) else {
            continue;
        };

        let Some(entity) = chain.entity() else {
            continue;
        };

        entities.insert(entity);
    }

    entities
}

/// Selects atoms belonging to one of `residues`.
///
/// Chains without an entity are intentionally skipped to preserve the behavior
/// of the original hierarchy expansion helper.
fn select_residues(
    structure: &Structure,
    universe: &AtomSelection,
    residues: &HashSet<ResidueIndex>,
) -> AtomSelection {
    let mut positions = Vec::new();

    for chain in structure.data().chains() {
        if chain.entity().is_none() {
            continue;
        }

        for residue in chain.residues() {
            if !residues.contains(&residue.index()) {
                continue;
            }

            append_residue_atoms(residue, universe, &mut positions);
        }
    }

    AtomSelection::from_sorted(positions)
}

/// Selects atoms belonging to one of `chains`.
///
/// Membership is tested once per chain rather than once per residue, reducing
/// hash lookups for long polymer chains.
fn select_chains(
    structure: &Structure,
    universe: &AtomSelection,
    chains: &HashSet<ChainIndex>,
) -> AtomSelection {
    let mut positions = Vec::new();

    for chain in structure.data().chains() {
        if chain.entity().is_none() || !chains.contains(&chain.index()) {
            continue;
        }

        for residue in chain.residues() {
            append_residue_atoms(residue, universe, &mut positions);
        }
    }

    AtomSelection::from_sorted(positions)
}

/// Selects atoms belonging to one of `entities`.
///
/// Entity membership is tested once per chain. Runtime is linear in the
/// traversed hierarchy plus universe-membership checks.
fn select_entities(
    structure: &Structure,
    universe: &AtomSelection,
    entities: &HashSet<EntityIndex>,
) -> AtomSelection {
    let mut positions = Vec::new();

    for chain in structure.data().chains() {
        let Some(entity) = chain.entity() else {
            continue;
        };

        if !entities.contains(&entity) {
            continue;
        }

        for residue in chain.residues() {
            append_residue_atoms(residue, universe, &mut positions);
        }
    }

    AtomSelection::from_sorted(positions)
}

/// Appends universe-visible atom indices from one residue in structure order.
///
/// The function performs no allocation beyond growth of `positions`; appended
/// indices remain sorted because hierarchy traversal is ordered.
#[inline]
fn append_residue_atoms(
    residue: ResidueRef<'_>,
    universe: &AtomSelection,
    positions: &mut Vec<u32>,
) {
    positions.extend(
        residue
            .atoms()
            .map(molframe_core::structure::AtomRef::index)
            .map(molframe_core::index::AtomIndex::get)
            .filter(|atom| universe.contains(*atom)),
    );
}
