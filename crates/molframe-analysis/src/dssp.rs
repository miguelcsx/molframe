//! Explicit semantic-role extraction for the shared native DSSP 4 kernel.
//!
//! Automatic enrichment and explicit analysis use identical hydrogen-bond and
//! pattern rules. A sorted spatial grid limits work to local C-alpha neighbours;
//! cross-chain bridges are allowed within one model.

pub use molframe_chem::DsspOptions;
use molframe_chem::{DsspBackbone, PolymerAtomRole, dssp_from_backbones};
use molframe_core::SecondaryStructure;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::{ResidueRef, Structure};

/// Why a secondary-structure assignment could not be evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum DsspError {
    /// Semantic polymer roles were not attached by an explicit profile.
    #[error("secondary structure requires explicit CCD polymer atom-role annotations")]
    MissingRoleAnnotation,
    /// A parameter is non-finite or outside its meaningful domain.
    #[error("secondary-structure options are invalid")]
    InvalidOptions,
    /// A backbone role that must be unique appears more than once.
    #[error("residue {residue:?} has multiple atoms for polymer role {role}")]
    AmbiguousRole {
        /// Residue containing the ambiguous role.
        residue: ResidueIndex,
        /// Stable integer role representation.
        role: i64,
    },
}

/// One residue and its secondary-structure state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SseRecord {
    /// The residue described.
    pub residue: ResidueIndex,
    /// Its assigned state.
    pub kind: SecondaryStructure,
}

define_soa_table! {
    /// Native columnar storage for secondary-structure assignments.
    pub struct SseTable for SseRecord {
        /// Residue indices.
        residue: ResidueIndex,
        /// Assigned secondary-structure states.
        kind: SecondaryStructure,
    }
}

impl SseTable {
    /// Returns computed DSSP rows as placement input without recomputation.
    pub fn placements(&self) -> impl Iterator<Item = (ResidueIndex, SecondaryStructure)> + '_ {
        self.residue()
            .iter()
            .copied()
            .zip(self.kind().iter().copied())
    }
}

/// Assigns modern native DSSP states using explicitly annotated polymer roles.
///
/// All complete heavy backbones are evaluable, including proline and chain starts
/// that cannot donate. Missing heavy roles remain `Unknown`. Rows are in residue
/// index order. Cross-chain sheets are computed globally within each model, while
/// local helix, turn, bend and PPII patterns require continuous peptide spans.
///
/// # Errors
///
/// Rejects invalid numerical options, missing semantic-role annotations, and
/// multiple atoms assigned to one unique backbone role.
pub fn secondary_structure(
    structure: &Structure,
    options: &DsspOptions,
) -> Result<SseTable, DsspError> {
    if !options.is_valid() {
        return Err(DsspError::InvalidOptions);
    }
    if !crate::chemistry::has_polymer_roles(structure) {
        return Err(DsspError::MissingRoleAnnotation);
    }
    let mut backbones = vec![DsspBackbone::default(); structure.residue_count()];
    for model in structure.data().models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                backbones[residue.index().as_usize()] = DsspBackbone {
                    model: model.index().get(),
                    chain: chain.index().get(),
                    proline: residue.name() == Some("PRO"),
                    ca: role_position(structure, residue, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?,
                    nitrogen: role_position(structure, residue, PolymerAtomRole::PROTEIN_NITROGEN)?,
                    carbon: role_position(
                        structure,
                        residue,
                        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                    )?,
                    oxygen: role_position(
                        structure,
                        residue,
                        PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                    )?,
                };
            }
        }
    }
    Ok(dssp_from_backbones(&backbones, options)
        .map_err(|_| DsspError::InvalidOptions)?
        .into_iter()
        .enumerate()
        .filter_map(|(i, kind)| {
            u32::try_from(i).ok().map(|i| SseRecord {
                residue: ResidueIndex::new(i),
                kind,
            })
        })
        .collect())
}

fn role_position(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
) -> Result<Option<[f32; 3]>, DsspError> {
    let mut matches = residue.atoms().filter(|atom| {
        crate::chemistry::polymer_role(structure, atom.index().get())
            .is_some_and(|role| role.intersects(required))
    });
    let first = matches
        .next()
        .and_then(molframe_core::structure::AtomRef::position);
    if matches.next().is_some() {
        Err(DsspError::AmbiguousRole {
            residue: residue.index(),
            role: required.code(),
        })
    } else {
        Ok(first)
    }
}

#[cfg(test)]
#[path = "dssp_tests.rs"]
mod tests;
