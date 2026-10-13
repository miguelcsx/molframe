//! Installing the dictionary's intra-component bonds between observed atoms.
//!
//! Observed atoms are ordered by name once so each bond endpoint is a binary
//! search instead of a scan of the residue.

use crate::Component;
use molframe_core::bond::{BondProvenance, BondRecord};
use molframe_core::optional::OptionalSymbol;
use molframe_core::structure::AtomRef;

pub(super) fn add_component_bonds(
    component: &Component,
    atoms: &[AtomRef<'_>],
    output: &mut Vec<BondRecord>,
) {
    // Scanning the residue's atoms once per bond endpoint costs
    // `O(bonds x atoms^2)` name comparisons. Ordering the observed atoms by
    // name once turns each endpoint into a binary search, and the sort is
    // stable so atoms sharing a name keep their original relative order.
    let mut by_name: Vec<(&str, AtomRef<'_>)> = atoms
        .iter()
        .filter_map(|atom| atom.name().map(|name| (name, *atom)))
        .collect();
    by_name.sort_by_key(|(name, _)| *name);

    for bond in component.bonds.iter() {
        for (_, atom_a) in named(&by_name, &bond.atom_a) {
            for (_, atom_b) in named(&by_name, &bond.atom_b) {
                if alt_compatible(*atom_a, *atom_b) {
                    output.push(BondRecord {
                        atom_a: atom_a.index(),
                        atom_b: atom_b.index(),
                        order: bond.order,
                        provenance: BondProvenance::ChemicalComponentDictionary,
                        // The dictionary states this bond exists inside a
                        // component; it is not a `struct_conn` classification of
                        // the assembly.
                        kind: OptionalSymbol::NONE,
                    });
                }
            }
        }
    }
}

/// The contiguous run of name-ordered atoms called `name`.
fn named<'a, 'b>(
    by_name: &'a [(&'a str, AtomRef<'b>)],
    name: &str,
) -> &'a [(&'a str, AtomRef<'b>)] {
    let start = by_name.partition_point(|(observed, _)| *observed < name);
    let end = by_name.partition_point(|(observed, _)| *observed <= name);
    match by_name.get(start..end) {
        Some(run) => run,
        None => &[],
    }
}

pub(super) fn alt_compatible(atom_a: AtomRef<'_>, atom_b: AtomRef<'_>) -> bool {
    match (atom_a.alt_id(), atom_b.alt_id()) {
        (Some(a), Some(b)) => a.is_blank() || b.is_blank() || a == b,
        _ => true,
    }
}
