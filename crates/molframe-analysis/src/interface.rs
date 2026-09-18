//! Residues that form the interface between two chains.
//!
//! An interface residue is one with any atom within the cutoff of an atom in the
//! other chain. The residues of both chains that meet this test are returned
//! together, sorted, so the result names the whole contact patch rather than one
//! side of it.
//!
//! The atom pairs come from the shared spatial search restricted to the two
//! chains, so the cost tracks the size of the smaller chain's neighbourhood.

use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, ChainRef, Structure};
use molframe_core::{ExecutionContext, index::ResidueIndex};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, StructureSpatial,
    reduce_pairs_within_unsorted,
};

/// Returns the residues at the interface between the two named chains.
///
/// A chain is matched by its label or author identifier. When either chain is
/// absent or the two never approach within the cutoff, the result is empty.
/// Residue indices are sorted and unique across both chains.
///
/// Runs in `O(chain atoms · local density)` time.
///
/// # Errors
///
/// Returns [`SpatialError`] for a non-finite or negative cutoff.
pub fn chain_interface(
    structure: &Structure,
    first: &str,
    second: &str,
    cutoff: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<Vec<ResidueIndex>, SpatialError> {
    let first_atoms = chain_atoms(structure, first);
    let second_atoms = chain_atoms(structure, second);
    let positions = structure.positions();
    let left = AtomSelection::from_sorted(first_atoms);
    let right = AtomSelection::from_sorted(second_atoms);
    let query = PairQuery {
        positions,
        left: &left,
        right: &right,
        cutoff,
        options: SpatialSearchOptions::with_backend(backend),
        periodic: None,
        context,
    };
    let parts = reduce_pairs_within_unsorted(&query, Vec::new, |residues: &mut Vec<u32>, pair| {
        push_interface_residues(structure, pair, residues);
    })?;

    let mut residues: Vec<u32> = Vec::new();
    for part in parts {
        residues.extend(part);
    }

    Ok(finish_interface_residues(residues))
}

/// Returns interface residues through a structure-bound spatial resolver.
///
/// A compiled plan uses this entry point so repeated interface/contact
/// operations can share the resolver's bounded cell/k-d index cache.  The
/// resolver remains borrowed from the same immutable structure snapshot, which
/// makes the lifetime and coordinate generation relationship explicit.
///
/// # Errors
///
/// Returns a diagnostic when the resolver rejects the workload or an atom index
/// cannot be mapped back to the structure hierarchy.
pub fn chain_interface_with_spatial(
    structure: &Structure,
    first: &str,
    second: &str,
    cutoff: f32,
    backend: SpatialBackend,
    spatial: &StructureSpatial<'_>,
) -> Result<Vec<ResidueIndex>, molframe_core::diagnostic::Diagnostic> {
    let first_atoms = chain_atoms(structure, first);
    let second_atoms = chain_atoms(structure, second);
    let left = AtomSelection::from_sorted(first_atoms);
    let right = AtomSelection::from_sorted(second_atoms);
    let pairs = spatial.pairs_with_backend(&left, &right, cutoff, backend)?;

    Ok(interface_residues(structure, pairs))
}

fn interface_residues(
    structure: &Structure,
    pairs: Vec<molframe_spatial::NeighborPair>,
) -> Vec<ResidueIndex> {
    let mut residues = Vec::new();
    for pair in pairs {
        push_interface_residues(structure, pair, &mut residues);
    }
    finish_interface_residues(residues)
}

/// Records both endpoints' residues for one contact pair.
///
/// Splitting the accumulation from the finish lets a streaming query feed the
/// same reduction without first materialising a pair vector.
fn push_interface_residues(
    structure: &Structure,
    pair: molframe_spatial::NeighborPair,
    residues: &mut Vec<u32>,
) {
    if let Some(residue) = residue_of(structure, pair.first) {
        residues.push(residue);
    }
    if let Some(residue) = residue_of(structure, pair.second) {
        residues.push(residue);
    }
}

/// Orders and deduplicates accumulated residues.
///
/// The sort makes the result independent of the order pairs arrived in, so a
/// streaming query and a materialising one agree exactly.
fn finish_interface_residues(mut residues: Vec<u32>) -> Vec<ResidueIndex> {
    residues.sort_unstable();
    residues.dedup();
    residues.into_iter().map(ResidueIndex::new).collect()
}

/// Collects the atom indices belonging to a chain matched by either identifier.
fn chain_atoms(structure: &Structure, name: &str) -> Vec<u32> {
    let mut atoms = Vec::new();
    for chain in structure.data().chains() {
        if !matches_name(chain, name) {
            continue;
        }
        for residue in chain.residues() {
            for atom in residue.atoms() {
                atoms.push(atom.index().get());
            }
        }
    }
    atoms.sort_unstable();
    atoms
}

fn matches_name(chain: ChainRef<'_>, name: &str) -> bool {
    chain.label() == Some(name) || chain.auth_label() == Some(name)
}

fn residue_of(structure: &Structure, atom: u32) -> Option<u32> {
    structure
        .data()
        .atom(molframe_core::index::AtomIndex::new(atom))
        .and_then(AtomRef::residue)
        .map(|residue| residue.index().get())
}

#[cfg(test)]
#[path = "interface_tests.rs"]
mod tests;
