//! Explicit semantic-role extraction for the shared native DSSP 4 kernel.
//!
//! Automatic enrichment and explicit analysis use identical hydrogen-bond and
//! pattern rules. A sorted spatial grid limits work to local C-alpha neighbours;
//! cross-chain bridges are allowed within one model.

pub use molframe_chem::DsspOptions;
use molframe_chem::{DsspBackbone, PolymerAtomRole, dssp_from_backbones};
use molframe_core::index::ResidueIndex;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{ResidueRef, Structure};
use molframe_core::{AnalysisPolicy, SecondaryStructure};

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
    secondary_structure_with_policy(structure, options, &AnalysisPolicy::default())
}

/// Like [`secondary_structure`], but selects atoms with the policy's
/// alternate-conformation rule before looking up backbone roles.
///
/// With the default policy one self-consistent conformer is used, so a residue
/// with altloc copies of its backbone atoms is evaluated rather than reported as
/// ambiguous. Under `KeepAll` the duplicates remain and are reported. When the
/// policy cannot be resolved for the structure (for example a ragged ensemble)
/// every atom is kept.
///
/// # Errors
///
/// Same as [`secondary_structure`].
pub fn secondary_structure_with_policy(
    structure: &Structure,
    options: &DsspOptions,
    policy: &AnalysisPolicy,
) -> Result<SseTable, DsspError> {
    if !options.is_valid() {
        return Err(DsspError::InvalidOptions);
    }
    if !crate::chemistry::has_polymer_roles(structure) {
        return Err(DsspError::MissingRoleAnnotation);
    }
    let selection = altloc_selection(structure, policy);
    let selection = selection.as_ref();
    let mut backbones = vec![DsspBackbone::default(); structure.residue_count()];
    for model in structure.data().models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                backbones[residue.index().as_usize()] = DsspBackbone {
                    model: model.index().get(),
                    chain: chain.index().get(),
                    proline: !can_donate(structure, residue, selection),
                    ca: role_position(
                        structure,
                        residue,
                        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                        selection,
                    )?,
                    nitrogen: role_position(
                        structure,
                        residue,
                        PolymerAtomRole::PROTEIN_NITROGEN,
                        selection,
                    )?,
                    carbon: role_position(
                        structure,
                        residue,
                        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                        selection,
                    )?,
                    oxygen: role_position(
                        structure,
                        residue,
                        PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                        selection,
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

/// The atoms kept by the policy's alternate-conformation rule, or `None` when
/// that rule cannot be resolved for this structure and every atom is kept.
pub(crate) fn altloc_selection(
    structure: &Structure,
    policy: &AnalysisPolicy,
) -> Option<AtomSelection> {
    structure.resolve_altlocs(policy).into_result().ok()
}

/// Whether the residue's backbone nitrogen can donate a hydrogen bond.
///
/// Decided from the structure, not the residue name. A backbone nitrogen with a
/// third heavy-atom neighbour besides the two chain atoms (the ring carbon of
/// proline, hydroxyproline and other proline-like residues) carries no N-H. The
/// neighbours come from the structure's bonds. The test counts heavy neighbours
/// inside the residue: an ordinary backbone nitrogen has only C-alpha there,
/// while a ring or N-alkylated nitrogen has two or more. The peptide bond to the
/// previous residue is deliberately not counted, because it exists only when
/// polymer links were built. Without bond data the same count is taken from
/// covalent-length contacts (under 1.9 angstrom) inside the residue, which is also
/// used when the bond table holds nothing for the nitrogen.
fn can_donate(
    structure: &Structure,
    residue: ResidueRef<'_>,
    selection: Option<&AtomSelection>,
) -> bool {
    let data = structure.data();
    let nitrogen = residue.atoms().find(|atom| {
        selection.is_none_or(|kept| kept.contains(atom.index().get()))
            && crate::chemistry::polymer_role(structure, atom.index().get())
                .is_some_and(|role| role.intersects(PolymerAtomRole::PROTEIN_NITROGEN))
    });
    let Some(nitrogen) = nitrogen else {
        return true;
    };
    let adjacency = data
        .bonds
        .is_available()
        .then(|| data.bonds.adjacency(structure.atom_count()));
    let bonded = adjacency
        .as_ref()
        .map_or(&[][..], |adjacency| adjacency.neighbours(nitrogen.index()));
    if bonded.is_empty() {
        return ring_contacts(residue, nitrogen, selection) < 2;
    }
    let heavy = bonded
        .iter()
        .filter(|neighbour| {
            data.atom(**neighbour).is_some_and(|atom| {
                atom.element() != Some(molframe_core::Element::HYDROGEN)
                    && atom
                        .residue()
                        .is_some_and(|owner| owner.index() == residue.index())
                    && selection.is_none_or(|kept| kept.contains(atom.index().get()))
            })
        })
        .count();
    heavy < 2
}

fn role_position(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
    selection: Option<&AtomSelection>,
) -> Result<Option<[f32; 3]>, DsspError> {
    let mut matches = residue.atoms().filter(|atom| {
        selection.is_none_or(|kept| kept.contains(atom.index().get()))
            && crate::chemistry::polymer_role(structure, atom.index().get())
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

/// Heavy atoms of the residue within covalent distance of its nitrogen.
fn ring_contacts(
    residue: ResidueRef<'_>,
    nitrogen: molframe_core::structure::AtomRef<'_>,
    selection: Option<&AtomSelection>,
) -> usize {
    const COVALENT_LIMIT: f32 = 1.9;
    let Some(origin) = nitrogen.position() else {
        return 0;
    };
    residue
        .atoms()
        .filter(|atom| {
            atom.index() != nitrogen.index()
                && atom.element() != Some(molframe_core::Element::HYDROGEN)
                && selection.is_none_or(|kept| kept.contains(atom.index().get()))
        })
        .filter_map(molframe_core::structure::AtomRef::position)
        .filter(|point| {
            let delta = [
                point[0] - origin[0],
                point[1] - origin[1],
                point[2] - origin[2],
            ];
            (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt()
                < COVALENT_LIMIT
        })
        .count()
}

#[cfg(test)]
#[path = "dssp_tests.rs"]
mod tests;
