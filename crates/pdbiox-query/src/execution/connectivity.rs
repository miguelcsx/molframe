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

/// Selects atoms within at most `depth` bonds of `target`, restricted to `universe`.
///
/// The bond graph is traversed breadth-first. Runtime is `O(V + E)` in the
/// visited subgraph and auxiliary space is `O(V + F)`, where the visited bitmap
/// uses one bit per atom and `F` is the maximum BFS frontier.
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

    if target.is_empty() {
        return Ok(AtomSelection::Empty);
    }

    let atom_count = structure.atom_count();
    let graph = structure.data().bonds.adjacency(atom_count);
    let atom_capacity = usize::try_from(atom_count).map_err(|_| {
        Diagnostic::new(Code::E4003).with_message("atom count exceeds platform capacity")
    })?;
    let mut seen = vec![false; atom_capacity];
    let mut queue = VecDeque::new();

    seed_queue(target, &mut seen, &mut queue);

    for _ in 0..depth {
        let breadth = queue.len();

        if breadth == 0 {
            break;
        }

        for _ in 0..breadth {
            let Some(atom) = queue.pop_front() else {
                break;
            };

            for neighbour in graph.neighbours(atom) {
                let Some(slot) = seen.get_mut(neighbour.as_usize()) else {
                    continue;
                };

                if *slot {
                    continue;
                }

                *slot = true;
                queue.push_back(*neighbour);
            }
        }
    }

    Ok(selection_from_seen(atom_count, universe, &seen))
}

/// Selects complete connected bond-graph components intersecting `target`.
///
/// Each reachable atom and bond edge is processed at most once, giving
/// `O(V + E)` runtime and `O(V + F)` auxiliary space.
pub(super) fn same_fragment(
    structure: &Structure,
    universe: &AtomSelection,
    target: &AtomSelection,
) -> Result<AtomSelection, Diagnostic> {
    require_graph(structure)?;

    if target.is_empty() {
        return Ok(AtomSelection::Empty);
    }

    let atom_count = structure.atom_count();
    let graph = structure.data().bonds.adjacency(atom_count);
    let atom_capacity = usize::try_from(atom_count).map_err(|_| {
        Diagnostic::new(Code::E4003).with_message("atom count exceeds platform capacity")
    })?;
    let mut seen = vec![false; atom_capacity];
    let mut queue = VecDeque::new();

    seed_queue(target, &mut seen, &mut queue);

    while let Some(atom) = queue.pop_front() {
        for neighbour in graph.neighbours(atom) {
            let Some(slot) = seen.get_mut(neighbour.as_usize()) else {
                continue;
            };

            if *slot {
                continue;
            }

            *slot = true;
            queue.push_back(*neighbour);
        }
    }

    Ok(selection_from_seen(atom_count, universe, &seen))
}

/// Seeds a BFS queue with valid atoms from `target`.
///
/// Duplicate or invalid atom indices do not enter the queue twice. Runtime is
/// `O(T)` for `T` target atoms and no allocation occurs beyond queue growth.
fn seed_queue(target: &AtomSelection, seen: &mut [bool], queue: &mut VecDeque<AtomIndex>) {
    for atom in target {
        let Some(slot) = seen.get_mut(AtomIndex::new(atom).as_usize()) else {
            continue;
        };

        if *slot {
            continue;
        }

        *slot = true;
        queue.push_back(AtomIndex::new(atom));
    }
}

/// Converts a dense visited bitmap into a sorted selection inside `universe`.
///
/// Runtime is `O(V * C_u)`, where `C_u` is the cost of
/// `AtomSelection::contains`, and output space is `O(M)` for `M` matches.
fn selection_from_seen(atom_count: u32, universe: &AtomSelection, seen: &[bool]) -> AtomSelection {
    let mut selected = Vec::new();

    for atom in 0..atom_count {
        if seen
            .get(AtomIndex::new(atom).as_usize())
            .copied()
            .is_some_and(|included| included)
            && universe.contains(atom)
        {
            selected.push(atom);
        }
    }

    AtomSelection::from_sorted(selected)
}

/// Verifies that the structure exposes a bond graph.
///
/// Returns `E4003` when connectivity information is unavailable. Runtime and
/// space are `O(1)`.
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
