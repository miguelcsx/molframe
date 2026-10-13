//! Donor and acceptor typing for structures that carry no hydrogens.
//!
//! Without a hydrogen a donor cannot be read from connectivity, so the role of
//! each heavy atom comes from what its residue is. That knowledge is a table, a
//! [`HeavyRoleTable`], and this module only applies it: one pass over the atoms,
//! asking the table for the role of each. An atom of a residue the table does not
//! class takes the generic rule: a nitrogen or oxygen can do both, unless the
//! structure carries donor and acceptor annotations, which are then authoritative.

use super::heavy_role_table::{ACCEPTOR, DONOR, HeavyRoleTable};
use molframe_core::column::Presence;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use molframe_core::{AtomAnnotation, AtomIndex};

const BOTH: u8 = DONOR | ACCEPTOR;

/// The donor and acceptor atoms of a structure.
pub(super) struct Roles {
    pub(super) donors: AtomSelection,
    pub(super) acceptors: AtomSelection,
}

/// Types every atom under `table`. `include_sulfur` admits the weak sulfur rows.
pub(super) fn assign_roles(
    structure: &Structure,
    table: &HeavyRoleTable,
    include_sulfur: bool,
) -> Roles {
    let count = structure.atom_count();
    let mut bits = vec![0u8; count as usize];
    let annotated = annotated_roles(structure);
    for chain in structure.data().chains() {
        let mut first = true;
        for residue in chain.residues() {
            let starts_chain = std::mem::replace(&mut first, false);
            for atom in residue.atoms() {
                let (Some(component), Some(name)) = (atom.component_name(), atom.name()) else {
                    continue;
                };
                let index = atom.index().get() as usize;
                let role = table
                    .role(component, name, starts_chain, include_sulfur)
                    .or_else(|| group_role(structure, atom.index(), annotated));
                if let (Some(role), Some(slot)) = (role, bits.get_mut(index)) {
                    *slot = role;
                }
            }
        }
    }
    let select = |mask: u8| {
        AtomSelection::from_sorted(
            (0..count)
                .filter(|&atom| bits.get(atom as usize).is_some_and(|bits| bits & mask != 0))
                .collect(),
        )
    };
    Roles {
        donors: select(DONOR),
        acceptors: select(ACCEPTOR),
    }
}

type RoleColumns<'a> = (
    &'a molframe_core::AnnotationColumn<bool>,
    &'a molframe_core::AnnotationColumn<bool>,
);

fn annotated_roles(structure: &Structure) -> Option<RoleColumns<'_>> {
    let annotations = structure.annotations();
    let (Some(AtomAnnotation::Boolean(donors)), Some(AtomAnnotation::Boolean(acceptors))) = (
        annotations.get(molframe_core::HBOND_DONOR_ANNOTATION),
        annotations.get(molframe_core::HBOND_ACCEPTOR_ANNOTATION),
    ) else {
        return None;
    };
    let count = structure.atom_count();
    (donors.len() == count && acceptors.len() == count).then_some((donors, acceptors))
}

/// Generic rule for atoms outside the standard residues: the annotations when
/// present, otherwise any nitrogen or oxygen can do both.
fn group_role(
    structure: &Structure,
    atom: AtomIndex,
    annotated: Option<RoleColumns<'_>>,
) -> Option<u8> {
    if let Some((donors, acceptors)) = annotated {
        let flag = |column: &molframe_core::AnnotationColumn<bool>, bit: u8| {
            if matches!(column.get(atom.get()), Some((true, Presence::Present))) {
                bit
            } else {
                0
            }
        };
        return Some(flag(donors, DONOR) | flag(acceptors, ACCEPTOR));
    }
    let element = structure.data().atom(atom)?.element()?;
    Some(if matches!(element.atomic_number(), 7 | 8) {
        BOTH
    } else {
        0
    })
}
