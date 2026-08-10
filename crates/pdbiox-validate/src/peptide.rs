//! Cis peptide bonds.
//!
//! The peptide bond is almost always trans: the two α-carbons lie on opposite
//! sides of the C–N bond and the ω torsion is near 180°. A cis bond, with ω
//! near 0° and the α-carbons on the same side, is rare and usually meaningful —
//! often a proline, sometimes an error — so it is worth flagging every one.
//!
//! ω is measured directly from CA(i)–C(i)–N(i+1)–CA(i+1) for consecutive
//! residues that are actually joined, judged by explicit connectivity when it
//! exists and otherwise by a C–N distance short enough to be a peptide bond.

use pdbiox_chem::PolymerAtomRole;
use pdbiox_core::Diagnostic;
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::{ResidueRef, Structure};

/// A residue whose peptide bond to the previous residue is cis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CisPeptide {
    /// The residue following the cis bond.
    pub residue: ResidueIndex,
    /// The ω torsion in degrees, near zero for a cis bond.
    pub omega: f64,
}

/// Finds every cis peptide bond, assigning it to the following residue.
///
/// A bond counts as cis when the magnitude of ω is at most `threshold_degrees`;
/// about 30° captures cis bonds without catching ordinary trans distortion.
/// Results are ordered by residue index.
///
/// Runs in `O(residues)` time.
///
/// # Errors
///
/// Returns a diagnostic when explicit polymer roles are absent or ambiguous.
pub fn cis_peptides(
    structure: &Structure,
    threshold_degrees: f64,
) -> Result<Vec<CisPeptide>, Diagnostic> {
    crate::backbone::require_polymer_roles(structure)?;
    let mut findings = Vec::new();
    for chain in structure.data().chains() {
        let residues: Vec<ResidueRef<'_>> = chain.residues().collect();
        for pair in residues.windows(2) {
            let current = pair[0];
            let next = pair[1];
            let Some(omega) = omega_between(structure, current, next)? else {
                continue;
            };
            let degrees = pdbiox_geom::degrees(omega);
            if degrees.abs() <= threshold_degrees {
                findings.push(CisPeptide {
                    residue: next.index(),
                    omega: degrees,
                });
            }
        }
    }
    findings.sort_by_key(|finding| finding.residue.get());
    Ok(findings)
}

/// Computes ω across the bond joining `current` to `next`, if they are joined.
fn omega_between(
    structure: &Structure,
    current: ResidueRef<'_>,
    next: ResidueRef<'_>,
) -> Result<Option<f64>, Diagnostic> {
    let alpha =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?;
    let carbon =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_CARBONYL_CARBON)?;
    let nitrogen = crate::backbone::role_atom(structure, next, PolymerAtomRole::PROTEIN_NITROGEN)?;
    let next_alpha =
        crate::backbone::role_atom(structure, next, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?;
    let (Some(alpha), Some(carbon), Some(nitrogen), Some(next_alpha)) =
        (alpha, carbon, nitrogen, next_alpha)
    else {
        return Ok(None);
    };
    if !crate::backbone::peptide_bonded(structure, carbon, nitrogen) {
        return Ok(None);
    }
    let (Some(alpha), Some(carbon), Some(nitrogen), Some(next_alpha)) = (
        alpha.position(),
        carbon.position(),
        nitrogen.position(),
        next_alpha.position(),
    ) else {
        return Ok(None);
    };
    Ok(pdbiox_geom::dihedral(alpha, carbon, nitrogen, next_alpha))
}

#[cfg(test)]
#[path = "peptide_tests.rs"]
mod tests;
