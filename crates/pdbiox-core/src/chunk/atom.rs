//! A contiguous run of atoms, stored column by column.
//!
//! Positions are not here. They live in one contiguous buffer per frame and a
//! chunk names the range of it that belongs to the chunk, so the whole
//! coordinate set can be handed out as a pointer while a kernel still reads a
//! chunk's worth at a time.

use super::parent::ParentMapping;
use super::stats::AtomChunkStats;
use crate::column::{EncodedColumn, Presence, ValidityMask};
use crate::coords::CoordinateBlock;
use crate::element::Element;
use crate::index::ResidueIndex;
use crate::optional::OptionalSymbol;
use crate::symbol::{AltId, SymbolId};
use crate::topology::ResidueTable;
use std::ops::Range;

/// One atom's worth of everything a chunk stores.
#[derive(Clone, Copy, Debug)]
pub struct AtomRecord {
    /// The atom's position, or `None` when the file did not record one.
    pub position: Option<[f32; 3]>,
    /// The element.
    pub element: Element,
    /// The normalised atom name.
    pub atom_name: SymbolId,
    /// The depositor's atom name, where it differs or where the file carried it.
    pub auth_atom_name: OptionalSymbol,
    /// The alternate-location label.
    pub alt_id: AltId,
    /// The residue this atom belongs to.
    pub residue: ResidueIndex,
    /// Occupancy, and whether it was recorded.
    pub occupancy: (f32, Presence),
    /// Temperature factor, and whether it was recorded.
    pub b_factor: (f32, Presence),
    /// Formal charge, and whether it was recorded.
    pub formal_charge: (i8, Presence),
    /// The file-local serial number, retained because other records cite it.
    pub atom_site_id: u32,
}

/// A contiguous run of atoms, stored column by column.
#[derive(Clone, Debug)]
pub struct AtomChunk {
    pub(super) atoms: Range<u32>,
    pub(super) model: u32,
    pub(super) element: EncodedColumn<u8>,
    pub(super) atom_name: EncodedColumn<SymbolId>,
    pub(super) auth_atom_name: Option<EncodedColumn<SymbolId>>,
    pub(super) alt_id: EncodedColumn<u32>,
    pub(super) residue: ParentMapping,
    pub(super) occupancy: EncodedColumn<f32>,
    pub(super) occupancy_validity: ValidityMask,
    pub(super) b_factor: EncodedColumn<f32>,
    pub(super) b_factor_validity: ValidityMask,
    pub(super) formal_charge: EncodedColumn<i8>,
    pub(super) formal_charge_validity: ValidityMask,
    pub(super) atom_site_id: EncodedColumn<u32>,
    pub(super) coord_validity: ValidityMask,
    pub(super) stats: AtomChunkStats,
}

impl AtomChunk {
    /// The number of atoms.
    #[must_use]
    pub fn len(&self) -> u32 {
        self.atoms.end.saturating_sub(self.atoms.start)
    }

    /// Returns true when the chunk holds no atoms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The atoms this chunk covers, in the structure's flat order.
    ///
    /// This is also the range of the coordinate buffer that belongs to it.
    #[must_use]
    pub fn atoms(&self) -> Range<u32> {
        self.atoms.clone()
    }

    /// What the chunk can say about itself without being read.
    #[must_use]
    pub const fn stats(&self) -> &AtomChunkStats {
        &self.stats
    }

    /// The model every atom in this chunk belongs to.
    #[must_use]
    pub const fn model(&self) -> u32 {
        self.model
    }

    /// The positions of this chunk's atoms, borrowed from the frame's buffer.
    #[must_use]
    pub fn positions<'a>(&self, coords: &'a CoordinateBlock) -> Option<&'a [[f32; 3]]> {
        coords.range(self.atoms.clone())
    }

    /// The element of the atom at `local`.
    #[must_use]
    pub fn element(&self, local: u32) -> Option<Element> {
        self.element.get(local).map(Element::from_atomic_number)
    }

    /// The normalised name of the atom at `local`.
    #[must_use]
    pub fn atom_name(&self, local: u32) -> Option<SymbolId> {
        self.atom_name.get(local)
    }

    /// The depositor's name for the atom at `local`, where the file carried one.
    #[must_use]
    pub fn auth_atom_name(&self, local: u32) -> Option<SymbolId> {
        self.auth_atom_name.as_ref()?.get(local)
    }

    /// The alternate-location label of the atom at `local`.
    #[must_use]
    pub fn alt_id(&self, local: u32) -> Option<AltId> {
        self.alt_id.get(local).map(AltId::from_raw)
    }

    /// The residue of the atom at `local`.
    #[must_use]
    pub fn residue(&self, local: u32, residues: &ResidueTable) -> Option<ResidueIndex> {
        self.residue.resolve(local, self.atoms.start, residues)
    }

    /// The occupancy of the atom at `local`, and whether it was recorded.
    #[must_use]
    pub fn occupancy(&self, local: u32) -> Option<(f32, Presence)> {
        Some((
            self.occupancy.get(local)?,
            self.occupancy_validity.get(local),
        ))
    }

    /// The temperature factor of the atom at `local`, and whether it was recorded.
    #[must_use]
    pub fn b_factor(&self, local: u32) -> Option<(f32, Presence)> {
        Some((self.b_factor.get(local)?, self.b_factor_validity.get(local)))
    }

    /// The formal charge of the atom at `local`, and whether it was recorded.
    #[must_use]
    pub fn formal_charge(&self, local: u32) -> Option<(i8, Presence)> {
        Some((
            self.formal_charge.get(local)?,
            self.formal_charge_validity.get(local),
        ))
    }

    /// The file-local serial number of the atom at `local`.
    #[must_use]
    pub fn atom_site_id(&self, local: u32) -> Option<u32> {
        self.atom_site_id.get(local)
    }

    /// Whether the position of the atom at `local` was recorded.
    #[must_use]
    pub fn has_position(&self, local: u32) -> bool {
        self.coord_validity.get(local).is_present()
    }
}
