//! Finding the residue an atom belongs to.
//!
//! The downward direction is free: a residue's atoms are a range, so it is two
//! loads and an addition. The upward direction is the awkward one, and it admits
//! a trade the other does not.
//!
//! Storing a residue number beside every atom answers it in one load and costs
//! four bytes per atom — sixteen megabytes on a four-million-atom virus capsid,
//! which is cache the geometric kernels do not get. Storing nothing answers it
//! by binary search over the residue offsets, which is fine occasionally and
//! poor in a loop.
//!
//! The default is neither: record the residue at the start of every block of a
//! hundred and twenty-eight atoms and scan forward within the block. That is
//! constant work per lookup at a thirty-second of the memory, because a block
//! spans a handful of residues and the scan never leaves one cache line.

use crate::index::ResidueIndex;
use crate::topology::ResidueTable;
use std::mem::size_of_val;

/// How an atom's residue is recovered.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub enum ParentMapping {
    /// Nothing stored; search the residue offsets.
    #[default]
    OffsetsOnly,
    /// One residue position per atom.
    Explicit(Vec<u32>),
    /// One residue position per block of atoms, refined by a bounded scan.
    BlockIndexed {
        /// Atoms per block.
        block: u16,
        /// The residue the first atom of each block belongs to.
        parents: Vec<u32>,
    },
}

impl ParentMapping {
    /// Atoms per block in the default mapping.
    ///
    /// A hundred and twenty-eight four-byte entries is one entry per five
    /// hundred and twelve bytes of atom data, and the scan it implies stays
    /// inside a few residues.
    pub const DEFAULT_BLOCK: u16 = 128;

    /// Builds the default mapping for a chunk whose atoms start at `first_atom`.
    ///
    /// `parent_of` is the residue of each atom in chunk order.
    #[must_use]
    pub fn block_indexed(parents_in_order: &[u32]) -> Self {
        let block = Self::DEFAULT_BLOCK;
        let parents = parents_in_order
            .iter()
            .step_by(usize::from(block))
            .copied()
            .collect();

        Self::BlockIndexed { block, parents }
    }

    /// Builds a mapping that answers in one load, for a chunk whose atoms are
    /// looked up individually often enough to justify the memory.
    #[must_use]
    pub fn explicit(parents_in_order: &[u32]) -> Self {
        Self::Explicit(parents_in_order.to_vec())
    }

    /// The bytes this mapping costs for `atoms` atoms.
    #[must_use]
    pub fn bytes(&self) -> usize {
        match self {
            Self::OffsetsOnly => 0,
            Self::Explicit(parents) | Self::BlockIndexed { parents, .. } => {
                size_of_val(parents.as_slice())
            }
        }
    }

    /// The residue of the atom at `local` within a chunk whose atoms begin at
    /// `first_atom` in the structure's flat order.
    ///
    /// `residues` is consulted only where the mapping cannot answer alone.
    #[must_use]
    pub fn resolve(
        &self,
        local: u32,
        first_atom: u32,
        residues: &ResidueTable,
    ) -> Option<ResidueIndex> {
        match self {
            Self::OffsetsOnly => Self::resolve_from_offsets(local, first_atom, residues),
            Self::Explicit(parents) => Self::resolve_explicit(local, parents),
            Self::BlockIndexed { block, parents } => {
                Self::resolve_block_indexed(local, first_atom, *block, parents, residues)
            }
        }
    }

    fn resolve_from_offsets(
        local: u32,
        first_atom: u32,
        residues: &ResidueTable,
    ) -> Option<ResidueIndex> {
        first_atom
            .checked_add(local)
            .and_then(|atom| residues.containing(atom))
    }

    fn resolve_explicit(local: u32, parents: &[u32]) -> Option<ResidueIndex> {
        usize::try_from(local)
            .ok()
            .and_then(|index| parents.get(index))
            .copied()
            .map(ResidueIndex::new)
    }

    fn resolve_block_indexed(
        local: u32,
        first_atom: u32,
        block: u16,
        parents: &[u32],
        residues: &ResidueTable,
    ) -> Option<ResidueIndex> {
        let start = Self::block_parent(local, block, parents)?;

        let atom = first_atom.checked_add(local)?;

        Self::scan_forward(start, atom, residues)
    }

    fn block_parent(local: u32, block: u16, parents: &[u32]) -> Option<u32> {
        let block = u32::from(block).max(1);
        let block_index = local / block;

        usize::try_from(block_index)
            .ok()
            .and_then(|index| parents.get(index))
            .copied()
    }

    /// Walks forward from a known residue until the one containing `atom`.
    ///
    /// Bounded by the number of residues a block spans, which is small by
    /// construction. Returns nothing rather than looping if the tables are
    /// inconsistent with each other.
    fn scan_forward(start: u32, atom: u32, residues: &ResidueTable) -> Option<ResidueIndex> {
        let start = usize::try_from(start).ok()?;

        (start..residues.len())
            .map_while(|candidate| {
                let candidate = u32::try_from(candidate).ok()?;

                let residue = ResidueIndex::new(candidate);
                let range = residues.atoms(residue)?;

                (atom >= range.start).then_some((residue, range))
            })
            .find_map(|(residue, range)| range.contains(&atom).then_some(residue))
    }
}

#[cfg(test)]
#[path = "parent_tests.rs"]
mod tests;
