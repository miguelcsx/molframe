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
    /// A component identity differing from the residue's primary identity.
    ///
    /// This occurs when alternate conformations model different residue types.
    /// Absence means the residue table supplies the component identity.
    pub alternate_component_id: OptionalSymbol,
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
    pub(super) auth_atom_name: Option<EncodedColumn<u32>>,
    pub(super) alternate_component_id: Option<EncodedColumn<u32>>,
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
    ///
    #[must_use]
    pub fn len(&self) -> u32 {
        self.atoms.end - self.atoms.start
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
        self.atom_range()
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
        coords.range(self.atom_range())
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

    /// Plain backing storage for atom names, when no decode is required.
    #[must_use]
    pub fn atom_names_plain(&self) -> Option<&[SymbolId]> {
        self.atom_name.as_slice()
    }

    /// The depositor's name for the atom at `local`, where the file carried one.
    #[must_use]
    pub fn auth_atom_name(&self, local: u32) -> Option<SymbolId> {
        optional_symbol(self.auth_atom_name.as_ref()?.get(local)?)
    }

    /// The atom-specific component identity, where it differs from its residue.
    #[must_use]
    pub fn alternate_component_id(&self, local: u32) -> Option<SymbolId> {
        optional_symbol(self.alternate_component_id.as_ref()?.get(local)?)
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
        with_presence(self.occupancy.get(local), &self.occupancy_validity, local)
    }

    /// Plain backing storage for occupancies, when no decode is required.
    #[must_use]
    pub fn occupancies_plain(&self) -> Option<&[f32]> {
        self.occupancy.as_slice()
    }

    /// The temperature factor of the atom at `local`, and whether it was recorded.
    #[must_use]
    pub fn b_factor(&self, local: u32) -> Option<(f32, Presence)> {
        with_presence(self.b_factor.get(local), &self.b_factor_validity, local)
    }

    /// Plain backing storage for temperature factors, when no decode is required.
    #[must_use]
    pub fn b_factors_plain(&self) -> Option<&[f32]> {
        self.b_factor.as_slice()
    }

    /// The formal charge of the atom at `local`, and whether it was recorded.
    #[must_use]
    pub fn formal_charge(&self, local: u32) -> Option<(i8, Presence)> {
        with_presence(
            self.formal_charge.get(local),
            &self.formal_charge_validity,
            local,
        )
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

    /// Reconstructs one logical row without changing its encoded storage.
    ///
    /// Editors use this only when a structural change requires rebuilding atom
    /// ordinals. Ordinary reads stay columnar and never materialise records.
    #[must_use]
    pub fn record(
        &self,
        local: u32,
        residues: &ResidueTable,
        position: Option<[f32; 3]>,
    ) -> Option<AtomRecord> {
        let position = if self.has_position(local) {
            Some(position?)
        } else {
            None
        };
        Some(AtomRecord {
            position,
            element: self.element(local)?,
            atom_name: self.atom_name(local)?,
            auth_atom_name: optional(self.auth_atom_name(local)),
            alternate_component_id: optional(self.alternate_component_id(local)),
            alt_id: self.alt_id(local)?,
            residue: self.residue(local, residues)?,
            occupancy: self.occupancy(local)?,
            b_factor: self.b_factor(local)?,
            formal_charge: self.formal_charge(local)?,
            atom_site_id: self.atom_site_id(local)?,
        })
    }

    fn atom_range(&self) -> Range<u32> {
        self.atoms.start..self.atoms.end
    }
}

fn with_presence<T>(
    value: Option<T>,
    validity: &ValidityMask,
    local: u32,
) -> Option<(T, Presence)> {
    value.map(|value| (value, validity.get(local)))
}

fn optional_symbol(raw: u32) -> Option<SymbolId> {
    (raw != u32::MAX).then(|| SymbolId::from_raw(raw))
}

fn optional(symbol: Option<SymbolId>) -> OptionalSymbol {
    match symbol {
        Some(symbol) => OptionalSymbol::some(symbol),
        None => OptionalSymbol::NONE,
    }
}
