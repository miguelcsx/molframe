//! Bond-graph selection operations.
//!
//! Traversal uses the core CSR cache. Each atom and edge is visited at most
//! once for a component expansion, while bounded bonded expansion stops after
//! the requested breadth-first depth.

use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::AtomIndex;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;

/// Number of bits stored in one visited-set machine word.
const WORD_BITS: usize = usize::BITS as usize;

/// Compact visited-atom set used by graph traversals.
///
/// Membership is stored as machine-word bits rather than one byte per atom,
/// reducing the working set while retaining direct atom-indexed access.
struct VisitSet {
    words: Vec<usize>,
    len: usize,
    marked: usize,
}

impl VisitSet {
    /// Creates an empty visited set capable of representing `len` atoms.
    fn new(len: usize) -> Self {
        let complete_words = len / WORD_BITS;
        let word_count = if len.is_multiple_of(WORD_BITS) {
            complete_words
        } else {
            complete_words + 1
        };

        Self {
            words: vec![0; word_count],
            len,
            marked: 0,
        }
    }

    /// Marks an atom and returns whether this is its first visit.
    ///
    /// Out-of-range indices are rejected without modifying the set.
    #[inline]
    fn mark_new(&mut self, index: usize) -> bool {
        if index >= self.len {
            return false;
        }

        let word_index = index / WORD_BITS;
        let bit_index = index % WORD_BITS;
        let mask = 1usize << bit_index;

        let Some(word) = self.words.get_mut(word_index) else {
            return false;
        };

        if *word & mask != 0 {
            return false;
        }

        *word |= mask;
        self.marked += 1;
        true
    }

    /// Returns whether an atom has already been visited.
    #[inline]
    fn contains(&self, index: usize) -> bool {
        if index >= self.len {
            return false;
        }

        let word_index = index / WORD_BITS;
        let bit_index = index % WORD_BITS;
        let mask = 1usize << bit_index;

        self.words
            .get(word_index)
            .is_some_and(|word| *word & mask != 0)
    }

    /// Returns the number of distinct atoms currently marked.
    #[inline]
    fn marked_count(&self) -> usize {
        self.marked
    }
}

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
    let (mut visited, mut queue) = traversal_state(atom_count, target)?;

    let mut head = 0usize;

    for _ in 0..depth {
        let level_end = queue.len();

        if head == level_end {
            break;
        }

        while head < level_end {
            let atom = queue[head];
            head += 1;

            enqueue_unseen(graph.neighbours(atom), &mut visited, &mut queue);
        }
    }

    Ok(selection_from_seen(universe, &visited))
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
    let (mut visited, mut queue) = traversal_state(atom_count, target)?;

    let mut head = 0usize;

    while head < queue.len() {
        let atom = queue[head];
        head += 1;

        enqueue_unseen(graph.neighbours(atom), &mut visited, &mut queue);
    }

    Ok(selection_from_seen(universe, &visited))
}

/// Creates the visited bitmap and BFS queue for one traversal.
///
/// Queue capacity is reserved from the target cardinality because every valid
/// seed must initially be retained until its adjacency has been processed.
fn traversal_state(
    atom_count: u32,
    target: &AtomSelection,
) -> Result<(VisitSet, Vec<AtomIndex>), Diagnostic> {
    let atom_capacity = usize::try_from(atom_count).map_err(|_| {
        Diagnostic::new(Code::E4003).with_message("atom count exceeds platform capacity")
    })?;

    let mut visited = VisitSet::new(atom_capacity);
    let queue_capacity = match usize::try_from(target.len()) {
        Ok(length) => length.min(atom_capacity),
        Err(_) => atom_capacity,
    };
    let mut queue = Vec::with_capacity(queue_capacity);

    seed_queue(target, &mut visited, &mut queue);

    Ok((visited, queue))
}

/// Seeds a BFS queue with valid atoms from `target`.
///
/// Duplicate or invalid atom indices do not enter the queue twice.
fn seed_queue(target: &AtomSelection, visited: &mut VisitSet, queue: &mut Vec<AtomIndex>) {
    for atom in target {
        let atom = AtomIndex::new(atom);

        if visited.mark_new(atom.as_usize()) {
            queue.push(atom);
        }
    }
}

/// Appends previously unseen graph neighbours to a traversal queue.
///
/// Invalid graph endpoints are ignored consistently with the existing
/// traversal behavior, while valid endpoints enter the queue only once.
#[inline]
fn enqueue_unseen<'a>(
    neighbours: impl IntoIterator<Item = &'a AtomIndex>,
    visited: &mut VisitSet,
    queue: &mut Vec<AtomIndex>,
) {
    for neighbour in neighbours {
        if visited.mark_new(neighbour.as_usize()) {
            queue.push(*neighbour);
        }
    }
}

/// Converts visited atoms into a sorted selection restricted to `universe`.
///
/// `AtomSelection` iteration already follows atom order, so the result remains
/// sorted without scanning unrelated structure atoms or performing a sort.
fn selection_from_seen(universe: &AtomSelection, visited: &VisitSet) -> AtomSelection {
    if universe.is_empty() || visited.marked_count() == 0 {
        return AtomSelection::Empty;
    }

    let seen = visited.marked_count();
    let selected_capacity = match usize::try_from(universe.len()) {
        Ok(length) => length.min(seen),
        Err(_) => seen,
    };
    let mut selected = Vec::with_capacity(selected_capacity);

    for atom in universe {
        if visited.contains(AtomIndex::new(atom).as_usize()) {
            selected.push(atom);
        }
    }

    if selected.is_empty() {
        AtomSelection::Empty
    } else {
        AtomSelection::from_sorted(selected)
    }
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
