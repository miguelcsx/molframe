//! Backbone connectivity shared by the peptide-bond and Ramachandran checks.
//!
//! Deciding whether two residues are actually joined is the same question
//! wherever a backbone torsion is measured, so it lives here once: explicit
//! connectivity is required. Callers whose source has no bonds must first apply
//! component chemistry under an explicit polymer-link policy.

use pdbiox_chem::PolymerAtomRole;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::{AtomAnnotation, Code, Diagnostic, Presence};

pub(crate) fn require_polymer_roles(structure: &Structure) -> Result<(), Diagnostic> {
    if matches!(
        structure
            .annotations()
            .get(pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION),
        Some(AtomAnnotation::Integer(_))
    ) {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4003)
            .with_context("required", "explicit CCD polymer atom-role annotations"))
    }
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

fn atom_role(structure: &Structure, atom: AtomRef<'_>) -> Option<PolymerAtomRole> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom.index().get())
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| PolymerAtomRole::from_code(code))
}

/// Whether a carbonyl carbon and the next residue's nitrogen form a bond.
pub(crate) fn peptide_bonded(
    structure: &Structure,
    carbon: AtomRef<'_>,
    nitrogen: AtomRef<'_>,
) -> bool {
    structure.data().bonds.is_available()
        && structure
            .data()
            .bonds
            .adjacency(structure.atom_count())
            .neighbours(carbon.index())
            .binary_search(&nitrogen.index())
            .is_ok()
}
