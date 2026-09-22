//! Small shared conversions used by hierarchy handles.

use crate::chunk::AtomChunk;
use crate::column::Presence;
use crate::symbol::SymbolId;
use std::ops::Range;

pub(super) fn range_or_empty(range: Option<Range<u32>>) -> Range<u32> {
    match range {
        Some(range) => range,
        None => 0..0,
    }
}

pub(super) fn recorded_value<T>(entry: (T, Presence)) -> Option<T> {
    let (value, presence) = entry;
    presence.is_present().then_some(value)
}

/// The first position in `range` whose atom carries `wanted`, if any.
///
/// One binary search locates the chunk holding the range's first atom; the walk
/// then advances at most once, because chunks tile the atom order in ascending
/// ranges and a chunk closes only before an atom that starts a new residue, so
/// a residue's own atoms never straddle a chunk boundary.
///
/// The search happens once per residue rather than once per atom. Repeating it
/// per atom is what a per-residue probe was paying for: every name read in a
/// three-atom residue re-searched the whole chunk list.
///
/// Returns `None` for a position outside every chunk, which a well-formed
/// structure cannot produce.
pub(super) fn find_named(chunks: &[AtomChunk], range: Range<u32>, wanted: SymbolId) -> Option<u32> {
    let mut cursor = chunks.partition_point(|chunk| chunk.atoms().end <= range.start);

    for position in range {
        loop {
            let chunk = chunks.get(cursor)?;
            let atoms = chunk.atoms();

            if position >= atoms.end {
                // The residue crossed into the next chunk, which the invariant
                // permits at most once for the whole range.
                cursor += 1;
                continue;
            }
            if position < atoms.start {
                return None;
            }
            if chunk.atom_name(position - atoms.start) == Some(wanted) {
                return Some(position);
            }
            break;
        }
    }

    None
}

#[cfg(test)]
#[path = "handle_helpers_tests.rs"]
mod tests;
