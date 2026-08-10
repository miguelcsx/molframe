//! Chain atoms and inter-set contacts, shared by the interface-based scores.
//!
//! Both `DockQ` and the quaternary-structure score start from the same two
//! questions: which atoms belong to a chain, and which pairs across two atom
//! sets are in contact. Answering them once here keeps the two scores agreeing
//! on what an interface contact is.

use pdbiox_core::contract::Namespace;
use pdbiox_core::structure::Structure;

use crate::CompareError;

/// The atom indices belonging to a chain in one explicit namespace, sorted.
pub(crate) fn chain_atoms(
    structure: &Structure,
    name: &str,
    namespace: Namespace,
) -> Result<Vec<usize>, CompareError> {
    let mut atoms = Vec::new();
    for chain in structure.data().chains() {
        let matches = match namespace {
            Namespace::Label => chain.label() == Some(name),
            Namespace::Auth => chain.auth_label() == Some(name),
            _ => return Err(CompareError::UnsupportedNamespace),
        };
        if !matches {
            continue;
        }
        for residue in chain.residues() {
            for atom in residue.atoms() {
                atoms.push(atom.index().as_usize());
            }
        }
    }
    atoms.sort_unstable();
    Ok(atoms)
}

/// The pairs `(a, b)` from the two index sets whose atoms are within `cutoff`.
pub(crate) fn contacts(
    positions: &[[f32; 3]],
    first: &[usize],
    second: &[usize],
    cutoff: f32,
) -> Vec<(usize, usize)> {
    let cutoff_squared = f64::from(cutoff) * f64::from(cutoff);
    let mut pairs = Vec::new();
    for &a in first {
        let Some(&point_a) = positions.get(a) else {
            continue;
        };
        for &b in second {
            let Some(&point_b) = positions.get(b) else {
                continue;
            };
            if pdbiox_geom::distance_squared(point_a, point_b) <= cutoff_squared {
                pairs.push((a, b));
            }
        }
    }
    pairs
}
