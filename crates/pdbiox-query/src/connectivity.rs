//! Bond-graph selection operations.
//!
//! Traversal uses the core CSR cache. Each atom and edge is visited at most
//! once for a component expansion, while bounded bonded expansion stops after
//! the requested breadth-first depth.

use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::index::AtomIndex;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use std::collections::VecDeque;

pub(super) fn bonded(
    structure: &Structure,
    universe: &AtomSelection,
    target: &AtomSelection,
    depth: u32,
) -> Result<AtomSelection, Diagnostic> {
    require_graph(structure)?;
    if depth == 0 {
        return Ok(universe.intersect(target));
    }
    let graph = structure.data().bonds.adjacency(structure.atom_count());
    let mut distance = vec![u32::MAX; structure.atom_count() as usize];
    let mut queue = VecDeque::new();
    for atom in target {
        if let Some(slot) = distance.get_mut(atom as usize) {
            *slot = 0;
            queue.push_back(AtomIndex::new(atom));
        }
    }
    while let Some(atom) = queue.pop_front() {
        let current = distance[atom.as_usize()];
        if current >= depth {
            continue;
        }
        for neighbour in graph.neighbours(atom) {
            let slot = &mut distance[neighbour.as_usize()];
            if *slot == u32::MAX {
                *slot = current + 1;
                queue.push_back(*neighbour);
            }
        }
    }
    let selected = distance
        .into_iter()
        .enumerate()
        .filter(|(atom, distance)| *distance <= depth && universe.contains(*atom as u32))
        .map(|(atom, _)| atom as u32)
        .collect();
    Ok(AtomSelection::from_sorted(selected))
}

pub(super) fn same_fragment(
    structure: &Structure,
    universe: &AtomSelection,
    target: &AtomSelection,
) -> Result<AtomSelection, Diagnostic> {
    require_graph(structure)?;
    let graph = structure.data().bonds.adjacency(structure.atom_count());
    let mut seen = vec![false; structure.atom_count() as usize];
    let mut queue = VecDeque::new();
    for atom in target {
        if let Some(slot) = seen.get_mut(atom as usize) {
            *slot = true;
            queue.push_back(AtomIndex::new(atom));
        }
    }
    while let Some(atom) = queue.pop_front() {
        for neighbour in graph.neighbours(atom) {
            if !seen[neighbour.as_usize()] {
                seen[neighbour.as_usize()] = true;
                queue.push_back(*neighbour);
            }
        }
    }
    let selected = seen
        .into_iter()
        .enumerate()
        .filter(|(atom, included)| *included && universe.contains(*atom as u32))
        .map(|(atom, _)| atom as u32)
        .collect();
    Ok(AtomSelection::from_sorted(selected))
}

fn require_graph(structure: &Structure) -> Result<(), Diagnostic> {
    if structure.data().bonds.is_available() {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4003).with_context("required", "bond graph"))
    }
}

#[cfg(test)]
#[path = "connectivity_tests.rs"]
mod tests;
