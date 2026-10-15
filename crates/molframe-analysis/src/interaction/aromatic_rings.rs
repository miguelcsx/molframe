//! Individual aromatic rings extracted from the bond graph.
//!
//! Aromatic atoms are annotated per atom, so a residue such as tryptophan or a
//! biaryl ligand carries several rings whose atoms are indistinguishable by
//! annotation alone. Fitting one plane to all of them describes none of the
//! rings. Here every ring is recovered as a cycle of the bond graph restricted
//! to aromatic atoms, and keeps its own atom set, centroid and normal.
//!
//! Rings are the smallest set of smallest rings: for each bond the shortest
//! cycle through it is found, and cycles are accepted in order of size when
//! they are independent of those already accepted (over the cycle space of the
//! graph). Fused systems therefore yield their constituent rings — indole gives
//! a five- and a six-membered ring — and never the envelope of the fusion.
//! Cycles longer than [`MAXIMUM_RING_SIZE`] are ignored, so the inner
//! macrocycle of a porphyrin does not appear as a ring of its own.
//!
//! A structure without a bond table has no ring topology, so it has no rings.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::structure::Structure;

/// Largest cycle, in atoms, that is treated as an aromatic ring.
pub const MAXIMUM_RING_SIZE: usize = 8;

/// One aromatic ring with its own plane.
#[derive(Clone, Debug, PartialEq)]
pub struct AromaticRing {
    /// The residue owning the ring's lowest-indexed atom.
    pub residue: ResidueIndex,
    /// The ring's atoms in ascending index order.
    pub atoms: Vec<AtomIndex>,
    /// The centroid of the ring atoms, in ångström.
    pub centroid: [f64; 3],
    /// The unit normal of the best-fit plane through the ring atoms.
    pub normal: [f64; 3],
}

impl AromaticRing {
    /// The lowest atom index of the ring, which identifies it within a structure.
    #[must_use]
    pub fn first_atom(&self) -> AtomIndex {
        match self.atoms.first() {
            Some(atom) => *atom,
            None => AtomIndex::new(0),
        }
    }
}

/// Extracts the aromatic rings of a structure, ordered by first atom.
///
/// Rings with any atom lacking a position, or fewer than three atoms, are
/// skipped. Runs in time linear in the number of aromatic bonds for bounded
/// ring size.
///
/// # Errors
///
/// Returns the eigen-solver error when a ring plane cannot be decomposed with
/// the selected numerical profile.
pub fn aromatic_rings(
    structure: &Structure,
    plane_fit: molframe_geom::EigenOptions,
) -> Result<Vec<AromaticRing>, molframe_geom::EigenError> {
    let mut rings = Vec::new();
    for atoms in ring_cycles(structure) {
        let points: Option<Vec<[f32; 3]>> = atoms
            .iter()
            .map(|atom| {
                structure
                    .data()
                    .atom(*atom)
                    .and_then(molframe_core::structure::AtomRef::position)
            })
            .collect();
        let Some(points) = points else {
            continue;
        };
        let Some(plane) = molframe_geom::best_fit_plane_with_options(&points, plane_fit)? else {
            continue;
        };
        let Some(residue) = structure
            .data()
            .atom(atoms[0])
            .and_then(molframe_core::structure::AtomRef::residue)
        else {
            continue;
        };
        rings.push(AromaticRing {
            residue: residue.index(),
            atoms,
            centroid: plane.centre,
            normal: plane.normal,
        });
    }
    Ok(rings)
}

/// Atom sets of the smallest set of smallest aromatic rings, ordered by first atom.
fn ring_cycles(structure: &Structure) -> Vec<Vec<AtomIndex>> {
    if !structure.data().bonds.is_available() {
        return Vec::new();
    }
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let aromatic = |atom: AtomIndex| crate::chemistry::is_aromatic(structure, atom.get());
    let neighbours = |atom: AtomIndex| -> Vec<AtomIndex> {
        adjacency
            .neighbours(atom)
            .iter()
            .copied()
            .filter(|other| aromatic(*other))
            .collect()
    };

    let mut edges: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    let mut atoms = BTreeSet::new();
    for atom in (0..structure.atom_count()).map(AtomIndex::new) {
        if !aromatic(atom) {
            continue;
        }
        for other in neighbours(atom) {
            if atom < other {
                let next = edges.len();
                edges.insert((atom.get(), other.get()), next);
                atoms.insert(atom.get());
                atoms.insert(other.get());
            }
        }
    }
    if edges.is_empty() {
        return Vec::new();
    }

    // Candidate cycles: the shortest cycle through each bond, as an edge set.
    let words = edges.len().div_ceil(64);
    let mut candidates: BTreeMap<(usize, Vec<u64>), Vec<AtomIndex>> = BTreeMap::new();
    for &(a, b) in edges.keys() {
        let Some(path) = shortest_path_avoiding(AtomIndex::new(a), AtomIndex::new(b), &neighbours)
        else {
            continue;
        };
        let mut vector = vec![0u64; words];
        for (position, atom) in path.iter().enumerate() {
            let next = path[(position + 1) % path.len()];
            let key = (atom.get().min(next.get()), atom.get().max(next.get()));
            if let Some(&index) = edges.get(&key) {
                vector[index / 64] ^= 1 << (index % 64);
            }
        }
        let mut members = path;
        members.sort_unstable();
        candidates.entry((members.len(), vector)).or_insert(members);
    }

    // Independent cycles in size order: rank of the cycle space is E - V + C.
    let rank = edges.len() + components(&atoms, &neighbours) - atoms.len();
    let mut basis: Vec<Vec<u64>> = Vec::new();
    let mut accepted = Vec::new();
    for ((_, vector), members) in candidates {
        if accepted.len() == rank {
            break;
        }
        if insert_independent(&mut basis, vector) {
            accepted.push(members);
        }
    }
    accepted.sort_unstable();
    accepted
}

/// Shortest path from `start` to `end` that does not use the direct bond, as atoms.
fn shortest_path_avoiding(
    start: AtomIndex,
    end: AtomIndex,
    neighbours: &impl Fn(AtomIndex) -> Vec<AtomIndex>,
) -> Option<Vec<AtomIndex>> {
    let mut parent: BTreeMap<AtomIndex, AtomIndex> = BTreeMap::new();
    let mut queue = VecDeque::from([(start, 1usize)]);
    parent.insert(start, start);
    while let Some((atom, size)) = queue.pop_front() {
        if size >= MAXIMUM_RING_SIZE {
            continue;
        }
        for next in neighbours(atom) {
            if atom == start && next == end {
                continue;
            }
            if parent.contains_key(&next) {
                continue;
            }
            parent.insert(next, atom);
            if next == end {
                let mut path = vec![end];
                let mut current = end;
                while current != start {
                    current = *parent.get(&current)?;
                    path.push(current);
                }
                return Some(path);
            }
            queue.push_back((next, size + 1));
        }
    }
    None
}

fn components(atoms: &BTreeSet<u32>, neighbours: &impl Fn(AtomIndex) -> Vec<AtomIndex>) -> usize {
    let mut seen = BTreeSet::new();
    let mut count = 0;
    for &atom in atoms {
        if !seen.insert(atom) {
            continue;
        }
        count += 1;
        let mut stack = vec![AtomIndex::new(atom)];
        while let Some(current) = stack.pop() {
            for next in neighbours(current) {
                if seen.insert(next.get()) {
                    stack.push(next);
                }
            }
        }
    }
    count
}

/// Adds `vector` to a GF(2) basis kept in echelon form; false when dependent.
fn insert_independent(basis: &mut Vec<Vec<u64>>, mut vector: Vec<u64>) -> bool {
    for row in basis.iter() {
        let Some(pivot) = lowest_bit(row) else {
            continue;
        };
        if vector[pivot / 64] >> (pivot % 64) & 1 == 1 {
            for (word, other) in vector.iter_mut().zip(row) {
                *word ^= other;
            }
        }
    }
    if lowest_bit(&vector).is_none() {
        return false;
    }
    basis.push(vector);
    true
}

fn lowest_bit(vector: &[u64]) -> Option<usize> {
    vector
        .iter()
        .enumerate()
        .find(|(_, word)| **word != 0)
        .map(|(index, word)| index * 64 + word.trailing_zeros() as usize)
}

#[cfg(test)]
#[path = "aromatic_rings_tests.rs"]
mod tests;
