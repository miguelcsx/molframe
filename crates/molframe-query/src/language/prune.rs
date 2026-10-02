//! Narrowing a scan with the statistics every chunk already carries.
//!
//! Each chunk records which elements it holds. A predicate that asks for one
//! element can therefore skip every chunk whose summary excludes it, without
//! reading a single atom of it. Only statistics that are exact are used:
//! `resid` is not pruned, because a chunk's residue bounds are hierarchy row
//! indices, not author sequence identifiers.

use molframe_core::Element;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;

/// The atoms of every chunk that may hold any of `elements`.
pub(crate) fn chunks_with_any_element(
    structure: &Structure,
    elements: &[Element],
) -> AtomSelection {
    AtomSelection::from_runs(
        structure
            .data()
            .chunks
            .iter()
            .filter(|chunk| {
                elements
                    .iter()
                    .any(|element| !chunk.stats().excludes_element(*element))
            })
            .map(molframe_core::AtomChunk::atoms),
    )
}

#[cfg(test)]
#[path = "prune_tests.rs"]
mod tests;
