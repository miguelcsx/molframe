//! Backbone connectivity shared by the peptide-bond and Ramachandran checks.
//!
//! Deciding whether two residues are actually joined is the same question
//! wherever a backbone torsion is measured, so it lives here once: explicit
//! connectivity decides when it exists, and a bounded carbon-nitrogen distance
//! is the fallback when the bond table does not list the link.

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
        .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom.index().get())
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| PolymerAtomRole::from_code(code))
}

/// Shortest carbonyl-carbon to amide-nitrogen distance accepted as a peptide
/// bond when the bond table is silent, in ångström. The ideal bond is about
/// 1.33 Å; the window leaves room for coordinate error without admitting a
/// non-bonded contact.
const PEPTIDE_BOND_MIN: f32 = 1.2;
/// Longest distance accepted as a peptide bond when the bond table is silent.
const PEPTIDE_BOND_MAX: f32 = 1.5;

/// Whether two atoms can belong to the same conformer: either has no alternate
/// location, or both carry the same one.
pub(crate) fn alt_compatible(a: AtomRef<'_>, b: AtomRef<'_>) -> bool {
    match (a.alt_id(), b.alt_id()) {
        (Some(a), Some(b)) => a.is_blank() || b.is_blank() || a == b,
        _ => true,
    }
}

/// Whether a carbonyl carbon and the next residue's nitrogen form a bond.
///
/// An entry in the bond table decides it. When the table has no such entry (or
/// there is no table) the pair still counts as joined if the two atoms are
/// altloc-compatible and lie within the peptide-bond distance window. Callers
/// pass consecutive residues of one chain, which is what keeps the distance
/// fallback from linking residues across chains.
pub(crate) fn peptide_bonded(
    structure: &Structure,
    carbon: AtomRef<'_>,
    nitrogen: AtomRef<'_>,
) -> bool {
    if structure.data().bonds.is_available()
        && structure
            .data()
            .bonds
            .adjacency(structure.atom_count())
            .neighbours(carbon.index())
            .binary_search(&nitrogen.index())
            .is_ok()
    {
        return true;
    }
    if !alt_compatible(carbon, nitrogen) {
        return false;
    }
    let (Some(c), Some(n)) = (carbon.position(), nitrogen.position()) else {
        return false;
    };
    let distance = c
        .iter()
        .zip(n)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f32>()
        .sqrt();
    (PEPTIDE_BOND_MIN..=PEPTIDE_BOND_MAX).contains(&distance)
}
