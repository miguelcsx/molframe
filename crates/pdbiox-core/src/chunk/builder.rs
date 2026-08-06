//! Building chunks from a stream of atoms in file order.
//!
//! A chunk is closed once it is at or past the target size *and* the next atom
//! starts a new residue. That is what keeps a chunk from straddling a residue
//! without splitting a large one, which in turn is what bounds the number of
//! chunks any residue-level operation has to touch.

use super::atom::{AtomChunk, AtomRecord};
use super::parent::ParentMapping;
use super::stats::AtomChunkStats;
use crate::column::{EncodedColumn, Presence, ValidityMask};
use crate::coords::CoordinateBlock;
use crate::symbol::SymbolId;

/// Atoms per chunk before a new one is started.
///
/// Chosen so the hot columns come to roughly a quarter of a megabyte, which is
/// what stays resident in second-level cache on the hardware this targets. It is
/// a tuned constant, not a promise, and nothing outside this crate depends on it.
pub const TARGET_CHUNK_ATOMS: u32 = 4096;

/// Accumulates atoms and closes chunks at residue boundaries.
///
/// Rows are appended in file order. A chunk is closed once it is at or past the
/// target size *and* the next atom starts a new residue, which is what keeps the
/// no-straddling guarantee without splitting a large residue.
///
/// # Examples
///
/// ```
/// use pdbiox_core::chunk::{AtomRecord, ChunkBuilder};
/// use pdbiox_core::{AltId, Element, OptionalSymbol, Presence, ResidueIndex, SymbolId};
///
/// let mut builder = ChunkBuilder::new();
/// builder.push(AtomRecord {
///     position: Some([1.0, 2.0, 3.0]),
///     element: Element::CARBON,
///     atom_name: SymbolId::from_raw(0),
///     auth_atom_name: OptionalSymbol::NONE,
///     alt_id: AltId::BLANK,
///     residue: ResidueIndex::new(0),
///     occupancy: (1.0, Presence::Present),
///     b_factor: (12.5, Presence::Present),
///     formal_charge: (0, Presence::Inapplicable),
///     atom_site_id: 1,
/// });
///
/// let (chunks, coords) = builder.finish();
/// assert_eq!(chunks.len(), 1);
/// assert_eq!(coords.len(), 1);
/// assert!(chunks[0].stats().elements.contains(Element::CARBON));
/// ```
#[derive(Debug, Default)]
pub struct ChunkBuilder {
    chunks: Vec<AtomChunk>,
    coords: CoordinateBlock,
    target: u32,
    model: u32,
    pending: Pending,
}

/// Column buffers for the chunk currently being filled.
#[derive(Debug, Default)]
struct Pending {
    first_atom: u32,
    element: Vec<u8>,
    atom_name: Vec<SymbolId>,
    auth_atom_name: Vec<SymbolId>,
    has_auth_names: bool,
    alt_id: Vec<u32>,
    residue: Vec<u32>,
    occupancy: Vec<f32>,
    occupancy_presence: Vec<Presence>,
    b_factor: Vec<f32>,
    b_factor_presence: Vec<Presence>,
    formal_charge: Vec<i8>,
    formal_charge_presence: Vec<Presence>,
    atom_site_id: Vec<u32>,
    coord_presence: Vec<Presence>,
    stats: AtomChunkStats,
}

impl Pending {
    fn len(&self) -> u32 {
        self.element.len() as u32
    }

    fn last_residue(&self) -> Option<u32> {
        self.residue.last().copied()
    }
}

impl ChunkBuilder {
    /// Creates a builder using the default chunk size.
    #[must_use]
    pub fn new() -> Self {
        Self {
            target: TARGET_CHUNK_ATOMS,
            ..Self::default()
        }
    }

    /// Creates a builder that closes chunks at `target` atoms.
    ///
    /// Exists for tests and for callers that know their access pattern; the
    /// default is what a file read uses.
    #[must_use]
    pub fn with_target(target: u32) -> Self {
        Self {
            target: target.max(1),
            ..Self::default()
        }
    }

    /// Reserves room for `atoms` positions.
    pub fn reserve(&mut self, atoms: usize) {
        self.coords = CoordinateBlock::with_capacity(atoms);
    }

    /// Declares which model the atoms that follow belong to.
    ///
    /// Closes the current chunk, because a chunk's model is constant by
    /// construction and that is what lets the column vanish entirely.
    pub fn start_model(&mut self, model: u32) {
        self.close();
        self.model = model;
    }

    /// Appends an atom.
    pub fn push(&mut self, record: AtomRecord) {
        let starts_new_residue = self.pending.last_residue() != Some(record.residue.get());
        if self.pending.len() >= self.target && starts_new_residue {
            self.close();
        }
        if self.pending.len() == 0 {
            self.pending.first_atom = self.coords.len();
        }

        if let Some(position) = record.position {
            self.coords.push(position);
            self.pending.stats.bounds.extend(position);
            self.pending.coord_presence.push(Presence::Present);
        } else {
            // The row still exists — an atom whose position was not recorded
            // is not an atom that was never modelled — so the buffer keeps a
            // slot and the validity mask carries the distinction.
            self.coords.push([f32::NAN; 3]);
            self.pending.coord_presence.push(Presence::Unknown);
            self.pending.stats.has_missing_coords = true;
        }

        self.pending.element.push(record.element.atomic_number());
        self.pending.atom_name.push(record.atom_name);
        match record.auth_atom_name.get() {
            Some(name) => {
                self.pending.auth_atom_name.push(name);
                self.pending.has_auth_names = true;
            }
            None => self.pending.auth_atom_name.push(record.atom_name),
        }
        self.pending.alt_id.push(record.alt_id.get());
        self.pending.residue.push(record.residue.get());
        self.pending.occupancy.push(record.occupancy.0);
        self.pending.occupancy_presence.push(record.occupancy.1);
        self.pending.b_factor.push(record.b_factor.0);
        self.pending.b_factor_presence.push(record.b_factor.1);
        self.pending.formal_charge.push(record.formal_charge.0);
        self.pending
            .formal_charge_presence
            .push(record.formal_charge.1);
        self.pending.atom_site_id.push(record.atom_site_id);

        self.observe(&record);
    }

    fn observe(&mut self, record: &AtomRecord) {
        let stats = &mut self.pending.stats;
        stats.elements.insert(record.element);
        stats.has_hydrogen |= record.element.is_hydrogen();
        stats.has_altloc |= !record.alt_id.is_blank();
        stats.model = self.model;
        if record.occupancy.1.is_present() {
            stats.occupancy.observe(record.occupancy.0);
        }
        if record.b_factor.1.is_present() {
            stats.b_factor.observe(record.b_factor.0);
        }
        let residue = record.residue.get();
        if self.pending.residue.len() == 1 {
            stats.residue_min = residue;
            stats.residue_max = residue;
        } else {
            stats.residue_min = stats.residue_min.min(residue);
            stats.residue_max = stats.residue_max.max(residue);
        }
    }

    /// Closes the chunk being filled, if it holds anything.
    fn close(&mut self) {
        if self.pending.len() == 0 {
            return;
        }
        let pending = std::mem::take(&mut self.pending);
        let len = pending.len();
        let first = pending.first_atom;

        let mut occupancy_validity: ValidityMask = pending.occupancy_presence.into_iter().collect();
        occupancy_validity.compact();
        let mut b_factor_validity: ValidityMask = pending.b_factor_presence.into_iter().collect();
        b_factor_validity.compact();
        let mut formal_charge_validity: ValidityMask =
            pending.formal_charge_presence.into_iter().collect();
        formal_charge_validity.compact();
        let mut coord_validity: ValidityMask = pending.coord_presence.into_iter().collect();
        coord_validity.compact();

        self.chunks.push(AtomChunk {
            atoms: first..first + len,
            model: self.model,
            element: EncodedColumn::encode(&pending.element),
            // Names are read by most selections, so they stay literal.
            atom_name: EncodedColumn::plain(&pending.atom_name),
            auth_atom_name: pending
                .has_auth_names
                .then(|| EncodedColumn::plain(&pending.auth_atom_name)),
            alt_id: EncodedColumn::encode(&pending.alt_id),
            residue: ParentMapping::block_indexed(&pending.residue),
            occupancy: EncodedColumn::plain(&pending.occupancy),
            occupancy_validity,
            b_factor: EncodedColumn::plain(&pending.b_factor),
            b_factor_validity,
            formal_charge: EncodedColumn::encode(&pending.formal_charge),
            formal_charge_validity,
            atom_site_id: EncodedColumn::encode(&pending.atom_site_id),
            coord_validity,
            stats: pending.stats,
        });
    }

    /// Closes the last chunk and returns everything built.
    #[must_use]
    pub fn finish(mut self) -> (Vec<AtomChunk>, CoordinateBlock) {
        self.close();
        (self.chunks, self.coords)
    }
}

#[cfg(test)]
#[path = "builder_tests.rs"]
mod tests;
