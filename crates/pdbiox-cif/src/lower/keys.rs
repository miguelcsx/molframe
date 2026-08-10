//! Deciding where one residue ends and the next begins.
//!
//! The file does not say. It gives every atom a set of annotations, and a
//! residue is a run of atoms that agree on them — so the whole question is which
//! annotations count.
//!
//! Two rules come from files that broke naive readers:
//!
//! - The component code is **excluded** from the key. One alternate location of
//!   a residue may carry a different component code, because a point mutation
//!   was modelled; that is one residue with two chemical identities, not two
//!   residues, and including the code in the key splits it.
//! - Where the annotations do not distinguish two adjacent groups at all, the
//!   caller must either reject the ambiguity or explicitly permit file-order
//!   inference. Real entries exist where consecutive residues carry identical
//!   annotations, and a reader that trusts the annotations merges them silently.

use pdbiox_core::optional::OptionalI32;

/// What identifies a residue, and therefore where its boundaries fall.
///
/// The component code is deliberately absent. Everything here is an interned
/// identifier or a number, so comparing two keys is integer work.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct ResidueKey {
    /// The model the atom belongs to.
    pub model: i64,
    /// The chain, interned.
    pub chain: u32,
    /// The sequence position, where the file gives one.
    pub label_seq: OptionalI32,
    /// The depositor's residue number.
    pub auth_seq: OptionalI32,
    /// The insertion code, interned, or absent.
    pub ins_code: u32,
}

impl ResidueKey {
    /// The value used for an interned field the file left empty.
    pub const ABSENT: u32 = u32::MAX;

    /// Returns true when this key says nothing that could distinguish it from
    /// `other` — meaning the annotations have run out and order must decide.
    #[must_use]
    pub fn is_indistinguishable_from(&self, other: &Self) -> bool {
        self == other
    }
}

/// Whether a new atom starts a new residue, and whether saying so needed a guess.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Boundary {
    /// The same residue continues.
    Same,
    /// A new residue begins, decided by the annotations.
    New,
    /// A new residue begins, decided by the order the file listed atoms in.
    ///
    /// Raised when the annotations cannot tell two adjacent groups apart. The
    /// atom names repeating is the signal: a residue does not contain two atoms
    /// of the same name in the same alternate location.
    NewByFileOrder,
}

/// Decides whether `next` continues the residue `current` describes.
///
/// `repeats_atom_name` says whether the incoming atom's name has already been
/// seen in the residue being filled, which is what betrays a boundary the
/// annotations do not mark.
#[must_use]
pub fn boundary(current: &ResidueKey, next: &ResidueKey, repeats_atom_name: bool) -> Boundary {
    if !current.is_indistinguishable_from(next) {
        return Boundary::New;
    }
    if repeats_atom_name {
        return Boundary::NewByFileOrder;
    }
    Boundary::Same
}

#[cfg(test)]
#[path = "keys_tests.rs"]
mod tests;
