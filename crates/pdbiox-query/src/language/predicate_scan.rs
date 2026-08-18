//! Ordered selection traversal shared by column predicates.

use super::AtomContext;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;

pub(super) struct SelectionCursor<I> {
    iter: I,
    next: Option<u32>,
}

impl<I> SelectionCursor<I>
where
    I: Iterator<Item = u32>,
{
    pub(super) fn new(mut iter: I) -> Self {
        let next = iter.next();
        Self { iter, next }
    }

    pub(super) fn is_exhausted(&self) -> bool {
        self.next.is_none()
    }

    pub(super) fn skip_before(&mut self, minimum: u32) {
        while self.next.is_some_and(|atom| atom < minimum) {
            self.next = self.iter.next();
        }
    }

    pub(super) fn matches(&mut self, atom: u32) -> bool {
        self.skip_before(atom);
        if self.next != Some(atom) {
            return false;
        }
        self.next = self.iter.next();
        true
    }

    pub(super) fn next_before(&mut self, end: u32) -> Option<u32> {
        let atom = self.next?;
        if atom >= end {
            return None;
        }
        self.next = self.iter.next();
        Some(atom)
    }
}

pub(crate) fn scan(
    structure: &Structure,
    universe: &AtomSelection,
    mut accepts: impl FnMut(AtomContext<'_>) -> bool,
) -> AtomSelection {
    if universe.is_empty() {
        return AtomSelection::Empty;
    }
    let mut selected = Vec::new();
    visit(structure, universe, |context| {
        if accepts(context) {
            selected.push(context.atom.index().get());
        }
    });
    selection_from_sorted(selected)
}

pub(super) fn visit<'a>(
    structure: &'a Structure,
    universe: &AtomSelection,
    mut visitor: impl FnMut(AtomContext<'a>),
) {
    let mut universe = SelectionCursor::new(universe.iter());
    if universe.is_exhausted() {
        return;
    }
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                let index = atom.index().get();
                if universe.matches(index) {
                    visitor(AtomContext {
                        atom,
                        residue,
                        chain,
                    });
                    if universe.is_exhausted() {
                        return;
                    }
                }
            }
        }
    }
}

pub(super) fn selection_from_sorted(selected: Vec<u32>) -> AtomSelection {
    if selected.is_empty() {
        AtomSelection::Empty
    } else {
        AtomSelection::from_sorted(selected)
    }
}
