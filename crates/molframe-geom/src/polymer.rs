//! Polymer geometry over explicitly connected backbone residues.
//!
//! A missing atom or chain break yields an absent torsion. It never fabricates
//! a peptide bond across a gap, and it never substitutes zero for an undefined
//! angle.

use crate::dihedral;

/// Backbone atoms needed for protein torsions.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BackboneResidue {
    /// Backbone nitrogen.
    pub nitrogen: Option<[f32; 3]>,
    /// Alpha carbon.
    pub alpha_carbon: Option<[f32; 3]>,
    /// Carbonyl carbon.
    pub carbonyl_carbon: Option<[f32; 3]>,
    /// Whether this residue has a peptide bond to the next input residue.
    pub connected_to_next: bool,
}

/// Protein backbone torsions for one residue, in radians.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BackboneTorsions {
    /// C(i-1)-N(i)-CA(i)-C(i).
    pub phi: Option<f64>,
    /// N(i)-CA(i)-C(i)-N(i+1).
    pub psi: Option<f64>,
    /// CA(i)-C(i)-N(i+1)-CA(i+1).
    pub omega: Option<f64>,
}

/// Computes φ, ψ and ω without crossing declared chain breaks.
///
/// Runs in `O(n)` time and returns one result per input residue.
#[must_use]
pub fn backbone_torsions(residues: &[BackboneResidue]) -> Vec<BackboneTorsions> {
    let mut torsions = Vec::with_capacity(residues.len());
    for (index, residue) in residues.iter().enumerate() {
        let previous = index
            .checked_sub(1)
            .and_then(|position| residues.get(position));
        let next = residues.get(index + 1);
        let phi = previous
            .filter(|previous| previous.connected_to_next)
            .and_then(|previous| {
                four(
                    previous.carbonyl_carbon,
                    residue.nitrogen,
                    residue.alpha_carbon,
                    residue.carbonyl_carbon,
                )
            });
        let (psi, omega) = if residue.connected_to_next {
            match next {
                Some(next) => (
                    four(
                        residue.nitrogen,
                        residue.alpha_carbon,
                        residue.carbonyl_carbon,
                        next.nitrogen,
                    ),
                    four(
                        residue.alpha_carbon,
                        residue.carbonyl_carbon,
                        next.nitrogen,
                        next.alpha_carbon,
                    ),
                ),
                None => (None, None),
            }
        } else {
            (None, None)
        };
        torsions.push(BackboneTorsions { phi, psi, omega });
    }
    torsions
}

fn four(
    first: Option<[f32; 3]>,
    second: Option<[f32; 3]>,
    third: Option<[f32; 3]>,
    fourth: Option<[f32; 3]>,
) -> Option<f64> {
    dihedral(first?, second?, third?, fourth?)
}

#[cfg(test)]
#[path = "polymer_tests.rs"]
mod tests;
