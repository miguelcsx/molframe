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
/// use molframe_core::chunk::{AtomRecord, ChunkBuilder};
/// use molframe_core::{AltId, Element, OptionalSymbol, Presence, ResidueIndex, SymbolId};
///
/// let mut builder = ChunkBuilder::new();
/// builder.push(AtomRecord {
///     position: Some([1.0, 2.0, 3.0]),
///     element: Element::CARBON,
///     atom_name: SymbolId::from_raw(0),
///     auth_atom_name: OptionalSymbol::NONE,
///     alternate_component_id: OptionalSymbol::NONE,
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
    auth_atom_name: Vec<u32>,
    has_auth_names: bool,
    alternate_component_id: Vec<u32>,
    has_alternate_components: bool,
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
    fn len(&self) -> usize {
        self.element.len()
    }

    fn is_empty(&self) -> bool {
        self.element.is_empty()
    }

    fn last_residue(&self) -> Option<u32> {
        self.residue.last().copied()
    }

    fn starts_new_residue(&self, residue: u32) -> bool {
        self.last_residue() != Some(residue)
    }

    fn should_close_before(&self, target: u32, residue: u32) -> bool {
        self.len() as u64 >= u64::from(target) && self.starts_new_residue(residue)
    }

    fn begin(&mut self, first_atom: u32, model: u32) {
        if !self.is_empty() {
            return;
        }

        self.first_atom = first_atom;
        self.stats.model = model;
    }

    fn push_record(&mut self, record: &AtomRecord) {
        self.element.push(record.element.atomic_number());
        self.atom_name.push(record.atom_name);
        self.push_auth_atom_name(record);
        self.push_alternate_component(record);
        self.alt_id.push(record.alt_id.get());
        self.residue.push(record.residue.get());

        push_with_presence(
            &mut self.occupancy,
            &mut self.occupancy_presence,
            record.occupancy,
        );

        push_with_presence(
            &mut self.b_factor,
            &mut self.b_factor_presence,
            record.b_factor,
        );

        push_with_presence(
            &mut self.formal_charge,
            &mut self.formal_charge_presence,
            record.formal_charge,
        );

        self.atom_site_id.push(record.atom_site_id);
        self.observe(record);
    }

    fn push_auth_atom_name(&mut self, record: &AtomRecord) {
        match record.auth_atom_name.get() {
            Some(name) => {
                self.auth_atom_name.push(name.get());
                self.has_auth_names = true;
            }
            None => self.auth_atom_name.push(u32::MAX),
        }
    }

    fn push_alternate_component(&mut self, record: &AtomRecord) {
        match record.alternate_component_id.get() {
            Some(component) => {
                self.alternate_component_id.push(component.get());
                self.has_alternate_components = true;
            }
            None => self.alternate_component_id.push(u32::MAX),
        }
    }

    fn observe(&mut self, record: &AtomRecord) {
        self.stats.elements.insert(record.element);
        self.stats.has_hydrogen |= record.element.is_hydrogen();
        self.stats.has_altloc |= !record.alt_id.is_blank();

        if record.occupancy.1.is_present() {
            self.stats.occupancy.observe(record.occupancy.0);
        }

        if record.b_factor.1.is_present() {
            self.stats.b_factor.observe(record.b_factor.0);
        }

        self.observe_residue(record.residue.get());
    }

    fn observe_residue(&mut self, residue: u32) {
        if self.residue.len() == 1 {
            self.stats.residue_min = residue;
            self.stats.residue_max = residue;
            return;
        }

        self.stats.residue_min = self.stats.residue_min.min(residue);
        self.stats.residue_max = self.stats.residue_max.max(residue);
    }

    fn reset(&mut self) {
        self.first_atom = 0;
        self.has_auth_names = false;
        self.has_alternate_components = false;

        self.element.clear();
        self.atom_name.clear();
        self.auth_atom_name.clear();
        self.alternate_component_id.clear();
        self.alt_id.clear();
        self.residue.clear();
        self.occupancy.clear();
        self.occupancy_presence.clear();
        self.b_factor.clear();
        self.b_factor_presence.clear();
        self.formal_charge.clear();
        self.formal_charge_presence.clear();
        self.atom_site_id.clear();
        self.coord_presence.clear();

        self.stats = AtomChunkStats::default();
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
        self.close_before(&record);
        self.pending.begin(self.coords.len(), self.model);
        self.push_position(record.position);
        self.pending.push_record(&record);
    }

    fn close_before(&mut self, record: &AtomRecord) {
        if self
            .pending
            .should_close_before(self.target, record.residue.get())
        {
            self.close();
        }
    }

    fn push_position(&mut self, position: Option<[f32; 3]>) {
        if let Some(position) = position {
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
    }

    /// Closes the chunk being filled, if it holds anything.
    fn close(&mut self) {
        if self.pending.is_empty() {
            return;
        }

        let chunk = self.build_chunk();
        self.chunks.push(chunk);
        self.pending.reset();
    }

    fn build_chunk(&mut self) -> AtomChunk {
        let end = self.coords.len();
        let pending = &mut self.pending;
        let first = pending.first_atom;
        let len = end - first;

        let occupancy_validity = compact_validity(&pending.occupancy_presence, len);
        let b_factor_validity = compact_validity(&pending.b_factor_presence, len);
        let formal_charge_validity = compact_validity(&pending.formal_charge_presence, len);
        let coord_validity = compact_validity(&pending.coord_presence, len);

        let element = EncodedColumn::encode(&pending.element);

        // Names are read by most selections, so they stay literal.
        let atom_name = EncodedColumn::plain(&pending.atom_name);

        let auth_atom_name = pending
            .has_auth_names
            .then(|| EncodedColumn::plain(&pending.auth_atom_name));
        let alternate_component_id = pending
            .has_alternate_components
            .then(|| EncodedColumn::encode(&pending.alternate_component_id));

        let alt_id = EncodedColumn::encode(&pending.alt_id);
        let residue = ParentMapping::block_indexed(&pending.residue);
        let occupancy = EncodedColumn::plain(&pending.occupancy);
        let b_factor = EncodedColumn::plain(&pending.b_factor);
        let formal_charge = EncodedColumn::encode(&pending.formal_charge);
        let atom_site_id = EncodedColumn::encode(&pending.atom_site_id);
        let stats = std::mem::take(&mut pending.stats);

        AtomChunk {
            atoms: first..end,
            model: self.model,
            element,
            atom_name,
            auth_atom_name,
            alternate_component_id,
            alt_id,
            residue,
            occupancy,
            occupancy_validity,
            b_factor,
            b_factor_validity,
            formal_charge,
            formal_charge_validity,
            atom_site_id,
            coord_validity,
            stats,
        }
    }

    /// Closes the last chunk and returns everything built.
    #[must_use]
    pub fn finish(mut self) -> (Vec<AtomChunk>, CoordinateBlock) {
        self.close();
        (self.chunks, self.coords)
    }
}

fn push_with_presence<T>(values: &mut Vec<T>, presences: &mut Vec<Presence>, entry: (T, Presence)) {
    let (value, presence) = entry;
    values.push(value);
    presences.push(presence);
}

fn compact_validity(presences: &[Presence], len: u32) -> ValidityMask {
    ValidityMask::from_bounded_iter(presences.iter().copied(), len)
}

#[cfg(test)]
#[path = "builder_tests.rs"]
mod tests;
