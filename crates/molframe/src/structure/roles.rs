//! Shared access to explicit polymer atom-role annotations.

use molframe_chem::PolymerAtomRole;
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use molframe_core::{AtomAnnotation, Code, Diagnostic, Presence};

pub(crate) fn require_polymer_roles(structure: &Structure) -> Result<(), Diagnostic> {
    if matches!(
        structure
            .annotations()
            .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION),
        Some(AtomAnnotation::Integer(_))
    ) {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4003)
            .with_context("required", "explicit CCD polymer atom-role annotations"))
    }
}

pub(crate) fn residue_has_role(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
) -> bool {
    residue
        .atoms()
        .any(|atom| atom_role(structure, atom).is_some_and(|role| role.intersects(required)))
}

pub(crate) fn role_atom<'a>(
    structure: &'a Structure,
    residue: ResidueRef<'a>,
    required: PolymerAtomRole,
) -> Result<Option<AtomRef<'a>>, Diagnostic> {
    let mut matching = residue
        .atoms()
        .filter(|atom| atom_role(structure, *atom).is_some_and(|role| role.intersects(required)));
    let first = matching.next();
    if matching.next().is_some() {
        Err(Diagnostic::new(Code::E4002)
            .with_context("residue", residue.index().to_string())
            .with_context(
                "required",
                format!("unique polymer role {}", required.code()),
            ))
    } else {
        Ok(first)
    }
}

pub(crate) fn atom_role(structure: &Structure, atom: AtomRef<'_>) -> Option<PolymerAtomRole> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom.index().get())
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| PolymerAtomRole::from_code(code))
}
