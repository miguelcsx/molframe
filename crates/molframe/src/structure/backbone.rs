//! Structure-aware projection into geometry kernels.

use molframe_chem::PolymerAtomRole;
use molframe_core::index::{AtomIndex, ChainIndex, ModelIndex, ResidueIndex};
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use molframe_core::{Diagnostic, Findings};
use molframe_geom::{BackboneResidue, BackboneTorsions};

use super::roles::{require_polymer_roles, residue_has_role, role_atom};

/// Backbone torsions associated with their structure residue.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BackboneTorsionRecord {
    /// Residue position in the structure topology.
    pub residue: ResidueIndex,
    /// Available φ, ψ and ω values in radians.
    pub torsions: BackboneTorsions,
}

/// One protein chain projected onto its CCD-annotated alpha-carbon trace.
#[derive(Clone, Debug, PartialEq)]
pub struct ProteinAlphaTrace {
    /// Chain position in the structure topology.
    pub chain: ChainIndex,
    /// One optional alpha-carbon position per protein-backbone residue.
    pub positions: Vec<Option<[f32; 3]>>,
}

/// Projects protein chains using CCD polymer atom roles rather than atom names.
///
/// # Errors
///
/// Returns a diagnostic when polymer-role annotations are absent or a residue
/// carries more than one alpha carbon.
pub fn structure_protein_alpha_traces(
    structure: &Structure,
) -> Result<Vec<ProteinAlphaTrace>, Findings> {
    require_polymer_roles(structure)?;
    let mut traces = Vec::new();
    for chain in structure.data().chains() {
        let mut positions = Vec::new();
        for residue in chain.residues().filter(|residue| {
            residue_has_role(structure, *residue, PolymerAtomRole::PROTEIN_BACKBONE)
        }) {
            positions.push(
                role_atom(structure, residue, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?
                    .and_then(AtomRef::position),
            );
        }
        if !positions.is_empty() {
            traces.push(ProteinAlphaTrace {
                chain: chain.index(),
                positions,
            });
        }
    }
    Ok(traces)
}

/// Computes protein backbone torsions chain by chain.
///
/// Explicit connectivity decides peptide continuity. Apply component chemistry
/// with a caller-selected polymer-link policy before this operation when the
/// source does not carry bonds.
///
/// # Errors
///
/// Returns a diagnostic when explicit polymer atom-role annotations are absent
/// or assign the same required semantic role to multiple atoms in one residue.
pub fn structure_backbone_torsions(
    structure: &Structure,
) -> Result<Vec<BackboneTorsionRecord>, Findings> {
    structure_backbone_torsions_model(structure, ModelIndex::new(0))
}

/// Computes protein backbone torsions for one dense model.
///
/// # Errors
///
/// Returns a diagnostic when explicit polymer atom-role annotations are absent
/// or ambiguous within a residue.
pub fn structure_backbone_torsions_model(
    structure: &Structure,
    model: ModelIndex,
) -> Result<Vec<BackboneTorsionRecord>, Findings> {
    require_polymer_roles(structure)?;
    let Some(positions) = structure.model_positions(model) else {
        return Ok(Vec::new());
    };
    let mut output = Vec::with_capacity(structure.residue_count());
    for chain in structure.data().chains() {
        let residues: Vec<ResidueRef<'_>> = chain
            .residues()
            .filter(|residue| {
                residue_has_role(structure, *residue, PolymerAtomRole::PROTEIN_BACKBONE)
            })
            .collect();
        let mut backbone = Vec::with_capacity(residues.len());
        for (position, residue) in residues.iter().enumerate() {
            let next = residues.get(position + 1).copied();
            backbone.push(project(structure, positions, *residue, next)?);
        }
        let torsions = molframe_geom::backbone_torsions(&backbone);
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
    Ok(output)
}

fn project(
    structure: &Structure,
    positions: &[[f32; 3]],
    residue: ResidueRef<'_>,
    next: Option<ResidueRef<'_>>,
) -> Result<BackboneResidue, Diagnostic> {
    let nitrogen = role_atom(structure, residue, PolymerAtomRole::PROTEIN_NITROGEN)?;
    let alpha_carbon = role_atom(structure, residue, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?;
    let carbonyl_carbon = role_atom(structure, residue, PolymerAtomRole::PROTEIN_CARBONYL_CARBON)?;
    Ok(BackboneResidue {
        nitrogen: nitrogen.and_then(|atom| position(positions, atom)),
        alpha_carbon: alpha_carbon.and_then(|atom| position(positions, atom)),
        carbonyl_carbon: carbonyl_carbon.and_then(|atom| position(positions, atom)),
        connected_to_next: next
            .map(|next| role_atom(structure, next, PolymerAtomRole::PROTEIN_NITROGEN))
            .transpose()?
            .flatten()
            .is_some_and(|next_nitrogen| {
                connected(structure, carbonyl_carbon, Some(next_nitrogen))
            }),
    })
}

fn connected(
    structure: &Structure,
    carbon: Option<AtomRef<'_>>,
    nitrogen: Option<AtomRef<'_>>,
) -> bool {
    let (Some(carbon), Some(nitrogen)) = (carbon, nitrogen) else {
        return false;
    };
    structure.data().bonds.is_available() && bonded(structure, carbon.index(), nitrogen.index())
}

fn position(positions: &[[f32; 3]], atom: AtomRef<'_>) -> Option<[f32; 3]> {
    positions.get(atom.index().as_usize()).copied()
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
#[path = "backbone_tests.rs"]
mod tests;
