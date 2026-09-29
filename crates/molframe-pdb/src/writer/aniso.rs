//! The anisotropic displacement record that follows an atom row.

use super::pdb::{alt_field, atom_name_field, insertion_code};
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use num_traits::ToPrimitive;
use std::fmt;

/// A tensor's identity columns, copied from the atom row before it.
pub(super) struct Record<'a> {
    /// The row the record pairs with, for its alternate-location field.
    pub atom: &'a AtomRef<'a>,
    /// The structure the atom belongs to, for resolving names.
    pub structure: &'a Structure,
    /// The residue the atom belongs to, for its sequence and insertion code.
    pub residue: &'a ResidueRef<'a>,
    /// The validated atom identifier the atom row carried.
    pub name: &'a str,
    /// The component identifier the atom row carried.
    pub component: &'a str,
    /// The chain label the atom row carried.
    pub chain: &'a str,
    /// The already-formatted sequence field of the atom row.
    pub seq: &'a str,
    /// The already-formatted serial field of the atom row.
    pub serial: &'a str,
    /// The tensor components, in `[U11, U22, U33, U12, U13, U23]` order.
    pub u: &'a [f32; 6],
}

impl Record<'_> {
    /// Writes the record under the identity columns it copied.
    pub(super) fn write(&self, out: &mut impl fmt::Write) {
        let units = self.u.map(field);
        let name = atom_name_field(self.name);
        let _ = writeln!(
            out,
            "ANISOU{serial}{name:<4}{alt:1}{comp:>3} {chain:>1}{seq}{ins:1} \
             {u0:>7}{u1:>7}{u2:>7}{u3:>7}{u4:>7}{u5:>7}  {element:>2}",
            serial = self.serial,
            alt = alt_field(self.structure, self.atom),
            comp = self.component,
            chain = self.chain,
            seq = self.seq,
            ins = insertion_code(*self.residue),
            element = match self.atom.element() {
                Some(element) => element.symbol(),
                None => "",
            },
            u0 = units[0],
            u1 = units[1],
            u2 = units[2],
            u3 = units[3],
            u4 = units[4],
            u5 = units[5],
        );
    }
}

/// One 7-column integer field, in units of 0.0001 ångström squared.
///
/// The reader divides by 10 000, so writing multiplies by it; a component the
/// format cannot hold as an integer clamps to the widest representable field
/// rather than silently becoming a different ellipsoid.
fn field(value: f32) -> i64 {
    let rounded = (f64::from(value) * 10_000.0).round();
    match rounded.to_i64() {
        Some(field) => field,
        None if rounded.is_sign_positive() => 99_999_999,
        None => -99_999_999,
    }
}
