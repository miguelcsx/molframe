//! Exact chemical-graph automorphism orbits.
//!
//! Colour refinement removes impossible correspondences in polynomial time;
//! the remaining exact backtracking is exponential in the worst case, as graph
//! automorphism itself can be. Molecular graphs are sparse and strongly
//! labelled, so refinement normally leaves only small symmetric groups.

use crate::numeric::usize_to_u32;
use crate::{Component, ComponentBond};
use pdbiox_core::contract::DictionaryVersion;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Atom-position orbits under element- and bond-labelled graph automorphisms.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EquivalenceClasses {
    classes: Arc<[Arc<[u32]>]>,
    class_by_atom: Arc<[u32]>,
}

impl EquivalenceClasses {
    /// Every orbit in order of its lowest atom position.
    #[must_use]
    pub fn classes(&self) -> &[Arc<[u32]>] {
        &self.classes
    }

    /// Orbit ordinal for a component atom position.
    #[must_use]
    pub fn class_of(&self, atom: u32) -> Option<u32> {
        self.class_by_atom.get(atom as usize).copied()
    }

    /// Whether two component atom positions are chemically interchangeable.
    #[must_use]
    pub fn equivalent(&self, atom_a: u32, atom_b: u32) -> bool {
        match (self.class_of(atom_a), self.class_of(atom_b)) {
            (Some(class_a), Some(class_b)) => class_a == class_b,
            _ => false,
        }
    }
}

/// Version-bound cache of equivalence classes by component identifier.
#[derive(Clone, Debug)]
pub struct EquivalenceCache {
    version: DictionaryVersion,
    entries: BTreeMap<Box<str>, Arc<EquivalenceClasses>>,
}

impl EquivalenceCache {
    /// Creates an empty cache for one exact CCD release.
    #[must_use]
    pub const fn new(version: DictionaryVersion) -> Self {
        Self {
            version,
            entries: BTreeMap::new(),
        }
    }

    /// Dictionary version whose component graphs key this cache.
    #[must_use]
    pub const fn version(&self) -> &DictionaryVersion {
        &self.version
    }

    /// Returns cached classes or computes and retains them once.
    pub fn get(&mut self, component: &Component) -> Arc<EquivalenceClasses> {
        if let Some(classes) = self.entries.get(component.id.as_ref()) {
            return classes.clone();
        }
        let classes = Arc::new(equivalence_classes(component));
        self.entries.insert(component.id.clone(), classes.clone());
        classes
    }

    /// Number of component graphs already analysed.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether no component graph has been analysed yet.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Computes exact automorphism orbits for a component graph.
#[must_use]
pub fn equivalence_classes(component: &Component) -> EquivalenceClasses {
    let graph = Graph::new(component);
    let count = graph.elements.len();
    let mut sets = DisjointSets::new(count);
    for source in 0..count {
        for target in source + 1..count {
            if graph.colours[source] == graph.colours[target]
                && graph.has_automorphism(source, target)
            {
                sets.union(source, target);
            }
        }
    }
    sets.finish()
}

/// Exact graph-automorphism enumeration exceeded its caller-selected bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AutomorphismLimit {
    /// Maximum complete mappings the caller allowed.
    pub limit: usize,
}

impl std::fmt::Display for AutomorphismLimit {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "component has more than {} graph automorphisms",
            self.limit
        )
    }
}

impl std::error::Error for AutomorphismLimit {}

/// Enumerates exact element- and bond-labelled component automorphisms.
///
/// Each mapping is indexed by source atom position and contains its target atom
/// position. Output is lexicographic and includes the identity mapping.
///
/// # Errors
///
/// Returns [`AutomorphismLimit`] instead of publishing a truncated set.
pub fn automorphisms(
    component: &Component,
    limit: usize,
) -> Result<Vec<Box<[u32]>>, AutomorphismLimit> {
    let graph = Graph::new(component);
    let mut forward = vec![usize::MAX; graph.elements.len()];
    let mut reverse = vec![usize::MAX; graph.elements.len()];
    let mut output = Vec::new();
    graph.enumerate(&mut forward, &mut reverse, 0, limit, &mut output)?;
    Ok(output)
}

#[derive(Debug)]
struct Graph {
    elements: Vec<u8>,
    edges: BTreeMap<(usize, usize), u8>,
    colours: Vec<u32>,
}

impl Graph {
    fn new(component: &Component) -> Self {
        let elements: Vec<u8> = component
            .atoms
            .iter()
            .map(|atom| atom.element.atomic_number())
            .collect();
        let mut edges = BTreeMap::new();
        let positions: BTreeMap<&str, usize> = component
            .atoms
            .iter()
            .enumerate()
            .map(|(position, atom)| (atom.name.as_ref(), position))
            .collect();
        for bond in component.bonds.iter() {
            let (Some(atom_a), Some(atom_b)) = (
                positions.get(bond.atom_a.as_ref()),
                positions.get(bond.atom_b.as_ref()),
            ) else {
                continue;
            };
            let label = bond_label(bond);
            edges.insert((*atom_a, *atom_b), label);
            edges.insert((*atom_b, *atom_a), label);
        }
        let colours = refine_colours(&elements, &edges);
        Self {
            elements,
            edges,
            colours,
        }
    }

    fn edge(&self, atom_a: usize, atom_b: usize) -> u8 {
        match self.edges.get(&(atom_a, atom_b)) {
            Some(edge) => *edge,
            None => 0,
        }
    }

    fn has_automorphism(&self, source: usize, target: usize) -> bool {
        let count = self.elements.len();
        let mut forward = vec![usize::MAX; count];
        let mut reverse = vec![usize::MAX; count];
        forward[source] = target;
        reverse[target] = source;
        self.extend(&mut forward, &mut reverse, 1)
    }

    fn extend(&self, forward: &mut [usize], reverse: &mut [usize], assigned: usize) -> bool {
        if assigned == forward.len() {
            return true;
        }
        let Some((source, candidates)) = self.next_candidates(forward, reverse) else {
            return false;
        };
        for target in candidates {
            forward[source] = target;
            reverse[target] = source;
            if self.extend(forward, reverse, assigned + 1) {
                return true;
            }
            forward[source] = usize::MAX;
            reverse[target] = usize::MAX;
        }
        false
    }

    fn enumerate(
        &self,
        forward: &mut [usize],
        reverse: &mut [usize],
        assigned: usize,
        limit: usize,
        output: &mut Vec<Box<[u32]>>,
    ) -> Result<(), AutomorphismLimit> {
        if assigned == forward.len() {
            if output.len() == limit {
                return Err(AutomorphismLimit { limit });
            }
            output.push(
                forward
                    .iter()
                    .map(|target| usize_to_u32(*target))
                    .collect::<Vec<_>>()
                    .into_boxed_slice(),
            );
            return Ok(());
        }
        let Some((source, candidates)) = self.next_candidates(forward, reverse) else {
            return Ok(());
        };
        for target in candidates {
            forward[source] = target;
            reverse[target] = source;
            self.enumerate(forward, reverse, assigned + 1, limit, output)?;
            forward[source] = usize::MAX;
            reverse[target] = usize::MAX;
        }
        Ok(())
    }

    fn next_candidates(&self, forward: &[usize], reverse: &[usize]) -> Option<(usize, Vec<usize>)> {
        let mut best: Option<(usize, Vec<usize>)> = None;
        for source in 0..forward.len() {
            if forward[source] != usize::MAX {
                continue;
            }
            let candidates: Vec<usize> = (0..reverse.len())
                .filter(|target| {
                    reverse[*target] == usize::MAX
                        && self.colours[source] == self.colours[*target]
                        && self.compatible(source, *target, forward)
                })
                .collect();
            if candidates.is_empty() {
                return None;
            }
            if best
                .as_ref()
                .is_none_or(|(_, current)| candidates.len() < current.len())
            {
                best = Some((source, candidates));
            }
        }
        best
    }

    fn compatible(&self, source: usize, target: usize, forward: &[usize]) -> bool {
        forward.iter().enumerate().all(|(other, mapped)| {
            *mapped == usize::MAX || self.edge(source, other) == self.edge(target, *mapped)
        })
    }
}

fn refine_colours(elements: &[u8], edges: &BTreeMap<(usize, usize), u8>) -> Vec<u32> {
    let count = elements.len();
    let mut colours: Vec<u32> = elements.iter().map(|element| u32::from(*element)).collect();
    loop {
        let signatures: Vec<(u32, Vec<(u8, u32)>)> = (0..count)
            .map(|atom| {
                let mut neighbours: Vec<(u8, u32)> = (0..count)
                    .filter_map(|other| {
                        let edge = match edges.get(&(atom, other)) {
                            Some(edge) => *edge,
                            None => 0,
                        };
                        (edge != 0).then_some((edge, colours[other]))
                    })
                    .collect();
                neighbours.sort_unstable();
                (colours[atom], neighbours)
            })
            .collect();
        let unique: BTreeSet<_> = signatures.iter().cloned().collect();
        let ids: BTreeMap<_, _> = unique
            .into_iter()
            .enumerate()
            .map(|(position, signature)| (signature, usize_to_u32(position)))
            .collect();
        let refined: Vec<u32> = signatures
            .iter()
            .filter_map(|signature| ids.get(signature).copied())
            .collect();
        if refined == colours {
            return colours;
        }
        colours = refined;
    }
}

const fn bond_label(bond: &ComponentBond) -> u8 {
    if bond.aromatic {
        return 5;
    }
    match bond.order {
        pdbiox_core::BondOrder::Unknown => 7,
        pdbiox_core::BondOrder::Single => 1,
        pdbiox_core::BondOrder::Double => 2,
        pdbiox_core::BondOrder::Triple => 3,
        pdbiox_core::BondOrder::Quadruple => 4,
        pdbiox_core::BondOrder::Aromatic => 5,
        pdbiox_core::BondOrder::Polymeric => 6,
    }
}

struct DisjointSets(Vec<usize>);

impl DisjointSets {
    fn new(count: usize) -> Self {
        Self((0..count).collect())
    }

    fn root(&self, mut atom: usize) -> usize {
        while self.0[atom] != atom {
            atom = self.0[atom];
        }
        atom
    }

    fn union(&mut self, atom_a: usize, atom_b: usize) {
        let (root_a, root_b) = (self.root(atom_a), self.root(atom_b));
        let (lower, upper) = if root_a <= root_b {
            (root_a, root_b)
        } else {
            (root_b, root_a)
        };
        self.0[upper] = lower;
    }

    fn finish(self) -> EquivalenceClasses {
        let mut grouped: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
        for atom in 0..self.0.len() {
            grouped
                .entry(self.root(atom))
                .or_default()
                .push(usize_to_u32(atom));
        }
        let classes: Vec<Arc<[u32]>> = grouped.into_values().map(Arc::from).collect();
        let mut class_by_atom = vec![0; self.0.len()];
        for (class, members) in classes.iter().enumerate() {
            for atom in members.iter().copied() {
                class_by_atom[atom as usize] = usize_to_u32(class);
            }
        }
        EquivalenceClasses {
            classes: classes.into(),
            class_by_atom: class_by_atom.into(),
        }
    }
}

#[cfg(test)]
#[path = "equivalence_tests.rs"]
mod tests;
