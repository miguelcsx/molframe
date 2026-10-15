//! Atom correspondence between two structures by chain, residue and atom name.
//!
//! Two entries of the same protein rarely number their atoms alike: a
//! different construct, a missing loop or an extra ligand shifts every later
//! index. What stays put is the name of the chain, the residue number and
//! name, and the atom name, so those identify the atom a comparison should
//! pair across the two.

use molframe_core::structure::Structure;
use std::collections::HashMap;

type AtomKey<'a> = (
    Option<&'a str>,
    Option<i32>,
    Option<&'a str>,
    Option<&'a str>,
    Option<&'a str>,
);

/// For each atom of `reference`, the atom of `target` with the same chain,
/// residue number, insertion code, residue name and atom name.
///
/// An atom with no such counterpart maps to `None`. When several atoms of one
/// structure share a key, as alternate conformations do, only the first takes
/// part and the rest map to `None`, so a correspondence never pairs two
/// reference atoms with one target atom.
#[must_use]
pub fn match_atoms_by_name(reference: &Structure, target: &Structure) -> Vec<Option<u32>> {
    let mut target_atoms: HashMap<AtomKey<'_>, u32> = HashMap::new();
    for (key, atom) in keyed_atoms(target) {
        target_atoms.entry(key).or_insert(atom);
    }

    let mut seen: std::collections::HashSet<AtomKey<'_>> = std::collections::HashSet::new();
    let mut matches = vec![None; reference.atom_count() as usize];
    for (key, atom) in keyed_atoms(reference) {
        if !seen.insert(key) {
            continue;
        }
        if let (Some(slot), Some(found)) = (matches.get_mut(atom as usize), target_atoms.get(&key))
        {
            *slot = Some(*found);
        }
    }
    matches
}

fn keyed_atoms(structure: &Structure) -> impl Iterator<Item = (AtomKey<'_>, u32)> {
    structure.data().chains().flat_map(|chain| {
        let chain_label = chain.auth_label().or_else(|| chain.label());
        chain.residues().flat_map(move |residue| {
            let number = residue.auth_seq_id().or_else(|| residue.label_seq_id());
            let insertion = residue.ins_code();
            let name = residue.name();
            residue.atoms().map(move |atom| {
                (
                    (chain_label, number, insertion, name, atom.name()),
                    atom.index().get(),
                )
            })
        })
    })
}

#[cfg(test)]
#[path = "native_match_tests.rs"]
mod tests;
