//! Ordered selection traversal shared by column predicates.

use super::AtomContext;
use molframe_core::chunk::AtomChunk;
use molframe_core::column::{Presence, ValidityMask};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{ChainRef, ResidueRef, Structure};

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
    // The universe's own count is the ceiling on what a predicate can accept,
    // so reserving it once replaces a growth sequence that reallocated as the
    // matches accumulated. A predicate accepting few atoms wastes the
    // reservation only until the buffer is dropped, while one accepting most of
    // the molecule — the common case for a broad selection — would otherwise
    // rebuild the buffer repeatedly as it grows.
    //
    // A universe larger than the address space reserves nothing and grows as it
    // did before: the reservation is an optimisation, and the selection it is
    // sized for cannot exist on a platform that cannot index it.
    let mut selected = match usize::try_from(universe.len()) {
        Ok(ceiling) => Vec::with_capacity(ceiling),
        Err(_) => Vec::new(),
    };
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

/// Selects the atoms of every residue `accepts`, deciding once per residue.
///
/// For a predicate that depends only on a residue and its chain this replaces a
/// decision per atom by one per residue, and returns the result as runs.
pub(crate) fn scan_residues<'a>(
    structure: &'a Structure,
    universe: &AtomSelection,
    mut accepts: impl FnMut(ChainRef<'a>, ResidueRef<'a>) -> bool,
) -> AtomSelection {
    if universe.is_empty() {
        return AtomSelection::Empty;
    }
    let residues = &structure.data().topology.residues;
    let mut runs = Vec::new();
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            if accepts(chain, residue)
                && let Some(atoms) = residues.atoms(residue.index())
            {
                runs.push(atoms);
            }
        }
    }
    let selected = AtomSelection::from_runs(runs);
    match universe {
        AtomSelection::All(_) => selected,
        _ => selected.intersect(universe),
    }
}

/// Selects the atoms `accepts`, reading each chunk's own columns directly.
///
/// An atom handle locates its chunk by binary search on every attribute read;
/// walking the chunks removes that search for predicates that read atom columns.
pub(crate) fn scan_chunks(
    structure: &Structure,
    universe: &AtomSelection,
    mut accepts: impl FnMut(&AtomChunk, u32) -> bool,
) -> AtomSelection {
    if universe.is_empty() {
        return AtomSelection::Empty;
    }
    let mut selected = match usize::try_from(universe.len()) {
        Ok(ceiling) => Vec::with_capacity(ceiling),
        Err(_) => Vec::new(),
    };
    let everything = matches!(universe, AtomSelection::All(_));
    let mut cursor = SelectionCursor::new(universe.iter());
    for chunk in structure.data().chunks.iter() {
        let atoms = chunk.atoms();
        for position in atoms.clone() {
            if !everything && !cursor.matches(position) {
                continue;
            }
            if accepts(chunk, position - atoms.start) {
                selected.push(position);
            }
        }
        if !everything && cursor.is_exhausted() {
            break;
        }
    }
    selection_from_sorted(selected)
}

/// How one `f32` atom column is read from a chunk.
pub(super) struct F32Column {
    /// The column's values as a plain slice, when no decode is needed.
    pub(super) plain: fn(&AtomChunk) -> Option<&[f32]>,
    /// Which positions of the column were recorded.
    pub(super) validity: fn(&AtomChunk) -> &ValidityMask,
    /// One value and whether it was recorded.
    pub(super) get: fn(&AtomChunk, u32) -> Option<(f32, Presence)>,
}

/// Selects the atoms whose recorded value of an `f32` column satisfies `test`.
///
/// A chunk that stores the column plainly with every value recorded is read as a
/// slice; any other chunk falls back to one lookup per atom.
pub(super) fn scan_f32(
    structure: &Structure,
    universe: &AtomSelection,
    column: &F32Column,
    test: impl Fn(f64) -> bool,
) -> AtomSelection {
    let everything = matches!(universe, AtomSelection::All(_));
    if !everything {
        return scan_chunks(structure, universe, |chunk, local| {
            (column.get)(chunk, local)
                .is_some_and(|(value, presence)| presence.is_present() && test(f64::from(value)))
        });
    }
    if universe.is_empty() {
        return AtomSelection::Empty;
    }
    let mut selected = Vec::new();
    for chunk in structure.data().chunks.iter() {
        let atoms = chunk.atoms();
        let plain = (column.plain)(chunk).filter(|values| {
            (column.validity)(chunk).is_all_present()
                && values.len() as u64 == u64::from(atoms.end - atoms.start)
        });
        if let Some(values) = plain {
            selected.extend(
                values
                    .iter()
                    .zip(atoms.clone())
                    .filter(|(value, _)| test(f64::from(**value)))
                    .map(|(_, position)| position),
            );
        } else {
            for position in atoms.clone() {
                let hit =
                    (column.get)(chunk, position - atoms.start).is_some_and(|(value, presence)| {
                        presence.is_present() && test(f64::from(value))
                    });
                if hit {
                    selected.push(position);
                }
            }
        }
    }
    selection_from_sorted(selected)
}
