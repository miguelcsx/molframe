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
            .step_by(block as usize)
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
                parents.len() * size_of::<u32>()
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
            Self::OffsetsOnly => residues.containing(first_atom.checked_add(local)?),
            Self::Explicit(parents) => parents.get(local as usize).copied().map(ResidueIndex::new),
            Self::BlockIndexed { block, parents } => {
                let block = u32::from(*block).max(1);
                let start = *parents.get((local / block) as usize)?;
                let atom = first_atom.checked_add(local)?;
                Self::scan_forward(start, atom, residues)
            }
        }
    }

    /// Walks forward from a known residue until the one containing `atom`.
    ///
    /// Bounded by the number of residues a block spans, which is small by
    /// construction. Returns nothing rather than looping if the tables are
    /// inconsistent with each other.
    fn scan_forward(start: u32, atom: u32, residues: &ResidueTable) -> Option<ResidueIndex> {
        let mut candidate = start;
        while (candidate as usize) < residues.len() {
            let range = residues.atoms(ResidueIndex::new(candidate))?;
            if range.contains(&atom) {
                return Some(ResidueIndex::new(candidate));
            }
            if atom < range.start {
                return None;
            }
            candidate += 1;
        }
        None
    }
}

#[cfg(test)]
#[path = "parent_tests.rs"]
mod tests;
