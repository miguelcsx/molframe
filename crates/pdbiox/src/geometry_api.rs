//! Structure-aware projection into geometry kernels.

use pdbiox_core::index::{AtomIndex, ResidueIndex};
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_geom::{BackboneResidue, BackboneTorsions};

const PEPTIDE_CUTOFF_SQUARED: f64 = 2.1 * 2.1;

/// Backbone torsions associated with their structure residue.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackboneTorsionRecord {
    /// Residue position in the structure topology.
    pub residue: ResidueIndex,
    /// Available φ, ψ and ω values in radians.
    pub torsions: BackboneTorsions,
}

/// Computes protein backbone torsions chain by chain.
///
/// Explicit connectivity decides peptide continuity when available. Otherwise
/// consecutive C-N atoms no further than 2.1 Å are treated as connected. The
/// geometry kernel remains structure-independent; this function owns only the
/// projection from hierarchy and bonds.
#[must_use]
pub fn structure_backbone_torsions(structure: &Structure) -> Vec<BackboneTorsionRecord> {
    let mut output = Vec::with_capacity(structure.residue_count());
    for chain in structure.data().chains() {
        let residues: Vec<ResidueRef<'_>> = chain.residues().collect();
        let mut backbone = Vec::with_capacity(residues.len());
        for (position, residue) in residues.iter().enumerate() {
            let next = residues.get(position + 1).copied();
            backbone.push(project(structure, *residue, next));
        }
        let torsions = pdbiox_geom::backbone_torsions(&backbone);
        output.extend(
            residues
                .into_iter()
                .zip(torsions)
                .map(|(residue, torsions)| BackboneTorsionRecord {
                    residue: residue.index(),
                    torsions,
                }),
        );
    }
    output
}

fn project(
    structure: &Structure,
    residue: ResidueRef<'_>,
    next: Option<ResidueRef<'_>>,
) -> BackboneResidue {
    let nitrogen = residue.atom("N");
    let alpha_carbon = residue.atom("CA");
    let carbonyl_carbon = residue.atom("C");
    BackboneResidue {
        nitrogen: nitrogen.and_then(AtomRef::position),
        alpha_carbon: alpha_carbon.and_then(AtomRef::position),
        carbonyl_carbon: carbonyl_carbon.and_then(AtomRef::position),
        connected_to_next: next
            .is_some_and(|next| connected(structure, carbonyl_carbon, next.atom("N"))),
    }
}

fn connected(
    structure: &Structure,
    carbon: Option<AtomRef<'_>>,
    nitrogen: Option<AtomRef<'_>>,
) -> bool {
    let (Some(carbon), Some(nitrogen)) = (carbon, nitrogen) else {
        return false;
    };
    if structure.data().bonds.is_available() {
        return bonded(structure, carbon.index(), nitrogen.index());
    }
    let (Some(carbon), Some(nitrogen)) = (carbon.position(), nitrogen.position()) else {
        return false;
    };
    pdbiox_geom::distance_squared(carbon, nitrogen) <= PEPTIDE_CUTOFF_SQUARED
}

fn bonded(structure: &Structure, carbon: AtomIndex, nitrogen: AtomIndex) -> bool {
    structure
        .data()
        .bonds
        .adjacency(structure.atom_count())
        .neighbours(carbon)
        .binary_search(&nitrogen)
        .is_ok()
}

#[cfg(test)]
#[path = "geometry_api_tests.rs"]
mod tests;
