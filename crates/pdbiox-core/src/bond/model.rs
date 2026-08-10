//! Compact connectivity storage and lazy graph traversal.
//!
//! Bonds live as sorted edge columns. Graph algorithms can request a cached
//! CSR view, built in O(atoms + bonds) time and storage, without imposing an
//! adjacency allocation on structures whose analyses never traverse bonds.

use crate::index::{AtomIndex, BondIndex};
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::{Arc, OnceLock};

/// The chemical order assigned to a bond.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondOrder {
    /// One shared electron pair.
    Single,
    /// Two shared electron pairs.
    Double,
    /// Three shared electron pairs.
    Triple,
    /// Four shared electron pairs.
    Quadruple,
    /// Delocalised aromatic connectivity.
    Aromatic,
    /// A linkage between consecutive polymer components.
    Polymeric,
    /// Connectivity is known but its order is not.
    #[default]
    Unknown,
}

/// Where a bond assignment came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondProvenance {
    /// An explicit record in the input structure.
    File,
    /// The chemical component dictionary.
    ChemicalComponentDictionary,
    /// An explicitly requested distance calculation.
    InferredDistance,
    /// An explicit caller edit.
    User,
}

/// One bond before it is packed into columns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BondRecord {
    /// First atom endpoint.
    pub atom_a: AtomIndex,
    /// Second atom endpoint.
    pub atom_b: AtomIndex,
    /// Chemical order.
    pub order: BondOrder,
    /// Source of this assignment.
    pub provenance: BondProvenance,
}

/// A compact, immutable edge table.
#[derive(Clone, Debug)]
pub struct BondTable {
    atom_a: Arc<Vec<AtomIndex>>,
    atom_b: Arc<Vec<AtomIndex>>,
    order: Arc<Vec<BondOrder>>,
    provenance: Arc<Vec<BondProvenance>>,
    adjacency: Arc<OnceLock<BondAdjacency>>,
    available: bool,
}

impl Default for BondTable {
    fn default() -> Self {
        Self {
            atom_a: Arc::new(Vec::new()),
            atom_b: Arc::new(Vec::new()),
            order: Arc::new(Vec::new()),
            provenance: Arc::new(Vec::new()),
            adjacency: Arc::new(OnceLock::new()),
            available: false,
        }
    }
}

impl BondTable {
    /// Number of unique undirected edges.
    #[must_use]
    pub fn len(&self) -> usize {
        self.atom_a.len()
    }

    /// Whether no connectivity is available.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atom_a.is_empty()
    }

    /// Whether connectivity was resolved, including a known empty graph.
    #[must_use]
    pub const fn is_available(&self) -> bool {
        self.available
    }

    /// Returns one edge by position.
    #[must_use]
    pub fn get(&self, bond: BondIndex) -> Option<BondRecord> {
        let position = bond.as_usize();
        Some(BondRecord {
            atom_a: *self.atom_a.get(position)?,
            atom_b: *self.atom_b.get(position)?,
            order: *self.order.get(position)?,
            provenance: *self.provenance.get(position)?,
        })
    }

    /// Walks records in stable endpoint order.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = BondRecord> + '_ {
        self.atom_a
            .iter()
            .zip(self.atom_b.iter())
            .zip(self.order.iter())
            .zip(self.provenance.iter())
            .map(|(((atom_a, atom_b), order), provenance)| BondRecord {
                atom_a: *atom_a,
                atom_b: *atom_b,
                order: *order,
                provenance: *provenance,
            })
    }

    /// Returns a cached compact adjacency index for graph traversal.
    #[must_use]
    pub fn adjacency(&self, atom_count: u32) -> &BondAdjacency {
        self.adjacency
            .get_or_init(|| BondAdjacency::build(self, atom_count))
    }
}

/// Deterministic bond-table construction with endpoint deduplication.
#[derive(Debug, Default)]
pub struct BondTableBuilder {
    records: BTreeMap<(u32, u32), BondRecord>,
}

impl BondTableBuilder {
    /// Creates an empty builder.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            records: BTreeMap::new(),
        }
    }

    /// Adds an undirected edge. The first assignment wins on duplicates.
    ///
    /// Self-edges are retained so structure validation can report them rather
    /// than silently changing the caller's data.
    pub fn push(&mut self, mut record: BondRecord) {
        if record.atom_b < record.atom_a {
            std::mem::swap(&mut record.atom_a, &mut record.atom_b);
        }
        self.records
            .entry((record.atom_a.get(), record.atom_b.get()))
            .or_insert(record);
    }

    /// Packs records into immutable columns.
    #[must_use]
    pub fn finish(self) -> BondTable {
        self.finish_with_availability(true)
    }

    /// Packs records while preserving whether the complete graph is known.
    ///
    /// A merge can retain explicit edges from one input while marking the
    /// combined graph unavailable because another input carried no connectivity.
    #[must_use]
    pub fn finish_with_availability(self, available: bool) -> BondTable {
        let mut atom_a = Vec::with_capacity(self.records.len());
        let mut atom_b = Vec::with_capacity(self.records.len());
        let mut order = Vec::with_capacity(self.records.len());
        let mut provenance = Vec::with_capacity(self.records.len());
        for record in self.records.into_values() {
            atom_a.push(record.atom_a);
            atom_b.push(record.atom_b);
            order.push(record.order);
            provenance.push(record.provenance);
        }
        BondTable {
            atom_a: Arc::new(atom_a),
            atom_b: Arc::new(atom_b),
            order: Arc::new(order),
            provenance: Arc::new(provenance),
            adjacency: Arc::new(OnceLock::new()),
            available,
        }
    }
}

/// A compressed sparse row view of an undirected bond graph.
#[derive(Clone, Debug, Default)]
pub struct BondAdjacency {
    offsets: Vec<u32>,
    neighbours: Vec<AtomIndex>,
}

impl BondAdjacency {
    fn build(table: &BondTable, atom_count: u32) -> Self {
        let mut offsets = vec![0u32; atom_count as usize + 1];
        for record in table.iter() {
            if record.atom_a.get() < atom_count && record.atom_b.get() < atom_count {
                offsets[record.atom_a.as_usize() + 1] += 1;
                offsets[record.atom_b.as_usize() + 1] += 1;
            }
        }
        for position in 1..offsets.len() {
            offsets[position] += offsets[position - 1];
        }
        let mut cursor = offsets.clone();
        let neighbour_count = match offsets.last() {
            Some(count) => *count,
            None => 0,
        };
        let mut neighbours = vec![AtomIndex::new(0); neighbour_count as usize];
        for record in table.iter() {
            if record.atom_a.get() >= atom_count || record.atom_b.get() >= atom_count {
                continue;
            }
            insert(&mut neighbours, &mut cursor, record.atom_a, record.atom_b);
            insert(&mut neighbours, &mut cursor, record.atom_b, record.atom_a);
        }
        for atom in 0..atom_count as usize {
            let range = offsets[atom] as usize..offsets[atom + 1] as usize;
            neighbours[range].sort_unstable();
        }
        Self {
            offsets,
            neighbours,
        }
    }

    /// Directly connected atoms in stable ordinal order.
    #[must_use]
    pub fn neighbours(&self, atom: AtomIndex) -> &[AtomIndex] {
        let Some(range) = self.range(atom) else {
            return &[];
        };
        &self.neighbours[range]
    }

    fn range(&self, atom: AtomIndex) -> Option<Range<usize>> {
        let next = atom.as_usize().checked_add(1)?;
        Some(*self.offsets.get(atom.as_usize())? as usize..*self.offsets.get(next)? as usize)
    }
}

fn insert(neighbours: &mut [AtomIndex], cursor: &mut [u32], atom: AtomIndex, value: AtomIndex) {
    let position = cursor[atom.as_usize()] as usize;
    if let Some(slot) = neighbours.get_mut(position) {
        *slot = value;
        cursor[atom.as_usize()] += 1;
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
