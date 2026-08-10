//! Shared access to CCD-derived atom annotations.

use pdbiox_core::annotation::AtomAnnotation;
use pdbiox_core::column::Presence;
use pdbiox_core::structure::Structure;

pub(crate) fn formal_charge(structure: &Structure, atom: u32) -> Option<i64> {
    let Some(AtomAnnotation::Integer(column)) = structure
        .annotations()
        .get(pdbiox_core::FORMAL_CHARGE_ANNOTATION)
    else {
        return None;
    };
    match column.get(atom) {
        Some((charge, Presence::Present)) => Some(charge),
        _ => None,
    }
}

pub(crate) fn component_kind(
    structure: &Structure,
    atom: u32,
) -> Option<pdbiox_chem::ComponentKind> {
    let Some(AtomAnnotation::Integer(column)) = structure
        .annotations()
        .get(pdbiox_core::COMPONENT_KIND_ANNOTATION)
    else {
        return None;
    };
    match column.get(atom) {
        Some((code, Presence::Present)) => pdbiox_chem::ComponentKind::from_code(code),
        _ => None,
    }
}

pub(crate) fn is_aromatic(structure: &Structure, atom: u32) -> bool {
    let Some(AtomAnnotation::Boolean(column)) = structure
        .annotations()
        .get(pdbiox_core::AROMATIC_ATOM_ANNOTATION)
    else {
        return false;
    };
    matches!(column.get(atom), Some((true, Presence::Present)))
}

pub(crate) fn has_polymer_roles(structure: &Structure) -> bool {
    matches!(
        structure
            .annotations()
            .get(pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION),
        Some(AtomAnnotation::Integer(_))
    )
}

pub(crate) fn polymer_role(
    structure: &Structure,
    atom: u32,
) -> Option<pdbiox_chem::PolymerAtomRole> {
    let Some(AtomAnnotation::Integer(column)) = structure
        .annotations()
        .get(pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION)
    else {
        return None;
    };
    match column.get(atom) {
        Some((code, Presence::Present)) => pdbiox_chem::PolymerAtomRole::from_code(code),
        _ => None,
    }
}
