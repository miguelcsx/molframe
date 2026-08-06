//! The residues of a structure.
//!
//! A residue is identified by its chain, its number, its insertion code and its
//! component. The insertion code is part of that identity rather than an
//! annotation on it: 163, 163A and 163B are three residues, and an antibody
//! numbering scheme depends on that entirely.

use crate::column::BitVec;
use crate::index::ResidueIndex;
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::symbol::SymbolId;
use std::ops::Range;

/// The residues of a structure.
#[derive(Clone, Debug, Default)]
pub struct ResidueTable {
    first_atom: Vec<u32>,
    atom_count: Vec<u32>,
    label_comp_id: Vec<SymbolId>,
    auth_comp_id: Vec<OptionalSymbol>,
    label_seq_id: Vec<OptionalI32>,
    auth_seq_id: Vec<OptionalI32>,
    ins_code: Vec<OptionalSymbol>,
    het: BitVec,
}

/// Everything a residue is created with.
#[derive(Clone, Copy, Debug)]
pub struct ResidueRecord {
    /// The normalised component code.
    pub label_comp_id: SymbolId,
    /// The depositor's component code, if the file carried one.
    pub auth_comp_id: OptionalSymbol,
    /// The sequence position within the entity. Absent for non-polymers.
    pub label_seq_id: OptionalI32,
    /// The depositor's residue number.
    pub auth_seq_id: OptionalI32,
    /// The insertion code, which is part of identity rather than an annotation.
    pub ins_code: OptionalSymbol,
    /// Whether the file recorded this residue as a heterogen.
    pub het: bool,
}

impl ResidueTable {
    /// The number of residues.
    #[must_use]
    pub fn len(&self) -> usize {
        self.label_comp_id.len()
    }

    /// Returns true when the structure has no residues.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.label_comp_id.is_empty()
    }

    /// Appends a residue covering a range of atoms.
    pub fn push(&mut self, record: ResidueRecord, atoms: Range<u32>) -> ResidueIndex {
        let position = self.label_comp_id.len() as u32;
        self.first_atom.push(atoms.start);
        self.atom_count.push(atoms.end.saturating_sub(atoms.start));
        self.label_comp_id.push(record.label_comp_id);
        self.auth_comp_id.push(record.auth_comp_id);
        self.label_seq_id.push(record.label_seq_id);
        self.auth_seq_id.push(record.auth_seq_id);
        self.ins_code.push(record.ins_code);
        self.het.push(record.het);
        ResidueIndex::new(position)
    }

    /// Replaces the atoms a residue covers.
    ///
    /// A reader that discovers a residue's extent only after reading its last
    /// atom needs this; nothing else should reach for it, because moving a
    /// residue's atoms without moving its neighbours' breaks the tiling.
    pub fn set_atoms(&mut self, residue: ResidueIndex, atoms: Range<u32>) {
        let Some(first) = self.first_atom.get_mut(residue.as_usize()) else {
            return;
        };
        *first = atoms.start;
        if let Some(count) = self.atom_count.get_mut(residue.as_usize()) {
            *count = atoms.end.saturating_sub(atoms.start);
        }
    }

    /// The atoms this residue contains.
    #[must_use]
    pub fn atoms(&self, residue: ResidueIndex) -> Option<Range<u32>> {
        let first = *self.first_atom.get(residue.as_usize())?;
        let count = *self.atom_count.get(residue.as_usize())?;
        Some(first..first.saturating_add(count))
    }

    /// The normalised component code.
    #[must_use]
    pub fn label_comp_id(&self, residue: ResidueIndex) -> Option<SymbolId> {
        self.label_comp_id.get(residue.as_usize()).copied()
    }

    /// The depositor's component code, if the file carried one.
    #[must_use]
    pub fn auth_comp_id(&self, residue: ResidueIndex) -> Option<SymbolId> {
        self.auth_comp_id
            .get(residue.as_usize())
            .and_then(|symbol| symbol.get())
    }

    /// The sequence position within the entity, absent for a non-polymer.
    #[must_use]
    pub fn label_seq_id(&self, residue: ResidueIndex) -> Option<i32> {
        self.label_seq_id
            .get(residue.as_usize())
            .and_then(|value| value.get())
    }

    /// The depositor's residue number.
    #[must_use]
    pub fn auth_seq_id(&self, residue: ResidueIndex) -> Option<i32> {
        self.auth_seq_id
            .get(residue.as_usize())
            .and_then(|value| value.get())
    }

    /// The insertion code.
    ///
    /// Part of identity: `163`, `163A` and `163B` are three residues, and an
    /// antibody numbering scheme depends on that entirely.
    #[must_use]
    pub fn ins_code(&self, residue: ResidueIndex) -> Option<SymbolId> {
        self.ins_code
            .get(residue.as_usize())
            .and_then(|symbol| symbol.get())
    }

    /// Whether the file recorded this residue as a heterogen.
    #[must_use]
    pub fn is_het(&self, residue: ResidueIndex) -> bool {
        self.het.test(residue.get())
    }

    /// The residue containing `atom`, found by searching the offsets.
    ///
    /// The offsets ascend, so this is a binary search. Chunks carry a denser
    /// mapping for the case where this question is asked per atom.
    #[must_use]
    pub fn containing(&self, atom: u32) -> Option<ResidueIndex> {
        let position = self.first_atom.partition_point(|first| *first <= atom);
        let candidate = ResidueIndex::new(position.checked_sub(1)? as u32);
        let range = self.atoms(candidate)?;
        range.contains(&atom).then_some(candidate)
    }
}
