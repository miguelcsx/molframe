//! Half-sphere exposure: how crowded each side of a residue is.
//!
//! A residue's α-carbon sits at the centre of a sphere split into two halves by
//! the plane through it perpendicular to the side-chain direction. The residues
//! whose α-carbons fall on the side-chain side are the upper count; those on the
//! other side are the lower count. The pair is a simple, rotation-free measure
//! of burial: a deeply buried residue is crowded on both sides, an exposed one
//! on neither.
//!
//! The side-chain direction is taken from the semantic alpha-carbon and
//! beta-carbon roles supplied by an explicit CCD-aware profile. A residue
//! without both roles is skipped rather than guessed at. Counts come from the
//! shared neighbour search over alpha carbons.

use pdbiox_chem::PolymerAtomRole;
use pdbiox_core::index::{AtomIndex, ResidueIndex};
use pdbiox_core::structure::{ResidueRef, Structure};
use pdbiox_core::{ExecutionContext, selection::AtomSelection};
use pdbiox_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, for_each_pairs_within_unsorted,
};

/// Why half-sphere exposure could not be evaluated.
#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
pub enum HseError {
    /// Semantic alpha/beta-carbon roles were not attached by an explicit profile.
    #[error("half-sphere exposure requires explicit CCD polymer atom-role annotations")]
    MissingRoleAnnotation,
    /// More than one atom in a residue carries a role that must be unique.
    #[error("residue {residue:?} has multiple atoms for polymer role {role}")]
    AmbiguousRole {
        /// Residue containing the ambiguous role.
        residue: ResidueIndex,
        /// Stable integer role representation.
        role: i64,
    },
    /// The selected spatial backend rejected the radius or coordinates.
    #[error("half-sphere neighbour search failed: {0}")]
    Spatial(#[from] SpatialError),
}

/// The crowding on each side of one residue.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HalfSphereExposure {
    /// The residue this describes.
    pub residue: ResidueIndex,
    /// Neighbouring α-carbons on the side-chain side.
    pub upper: u32,
    /// Neighbouring α-carbons on the opposite side.
    pub lower: u32,
}

/// One residue's α-carbon, its side-chain direction, and its running counts.
struct Centre {
    residue: ResidueIndex,
    atom: AtomIndex,
    alpha: [f32; 3],
    direction: [f32; 3],
    upper: u32,
    lower: u32,
}

/// Computes half-sphere exposure for residues whose explicit chemistry profile
/// identifies one alpha carbon and one beta carbon.
///
/// Two α-carbons count as neighbours when within `radius`, about 13 Å by
/// convention. Results are ordered by residue index.
///
/// Runs in `O(residues · local density)` time.
///
/// # Errors
///
/// Returns [`HseError::MissingRoleAnnotation`] when semantic atom roles were not
/// applied, [`HseError::AmbiguousRole`] for non-unique roles, and
/// [`HseError::Spatial`] for an invalid radius or spatial input.
pub fn half_sphere_exposure(
    structure: &Structure,
    radius: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<Vec<HalfSphereExposure>, HseError> {
    require_role_annotation(structure)?;
    let mut centres = Vec::new();
    let mut centre_of_atom = vec![None; structure.atom_count() as usize];
    for residue in structure.data().residues() {
        let (Some(alpha_index), Some(beta_index)) = (
            role_atom(structure, residue, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?,
            role_atom(structure, residue, PolymerAtomRole::PROTEIN_BETA_CARBON)?,
        ) else {
            continue;
        };
        let (Some(alpha), Some(beta)) = (
            structure.positions().get(alpha_index.as_usize()).copied(),
            structure.positions().get(beta_index.as_usize()).copied(),
        ) else {
            continue;
        };
        let direction = [beta[0] - alpha[0], beta[1] - alpha[1], beta[2] - alpha[2]];
        if let Some(slot) = centre_of_atom.get_mut(alpha_index.as_usize()) {
            *slot = Some(centres.len());
        }
        centres.push(Centre {
            residue: residue.index(),
            atom: alpha_index,
            alpha,
            direction,
            upper: 0,
            lower: 0,
        });
    }

    let alpha_atoms: Vec<u32> = centres.iter().map(|centre| centre.atom.get()).collect();
    let mut sorted = alpha_atoms;
    sorted.sort_unstable();
    let selection = AtomSelection::from_sorted(sorted);
    // Every pair collapses into two per-residue counters, so nothing is gained
    // by holding the pairs: they are classified as they are produced.
    //
    // Deliberately serial. A blocked reduction allocates one accumulator per
    // block, and this accumulator is one counter pair per residue — so parallel
    // execution would cost residues x blocks, which for a large structure is
    // far more memory than the search saves in time. A sparse per-block map
    // would fix the memory and add a hash lookup to every pair, which is most
    // of the per-pair work here. The decomposition this kernel wants is over
    // residues, not over cells.
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions: structure.positions(),
            left: &selection,
            right: &selection,
            cutoff: radius,
            options: SpatialSearchOptions::with_backend(backend),
            periodic: None,
            context,
        },
        |pair| {
            let (Some(&Some(a)), Some(&Some(b))) = (
                centre_of_atom.get(pair.first as usize),
                centre_of_atom.get(pair.second as usize),
            ) else {
                return;
            };
            classify(&mut centres, a, b);
            classify(&mut centres, b, a);
        },
    )?;

    let mut result: Vec<HalfSphereExposure> = centres
        .iter()
        .map(|centre| HalfSphereExposure {
            residue: centre.residue,
            upper: centre.upper,
            lower: centre.lower,
        })
        .collect();
    result.sort_by_key(|exposure| exposure.residue.get());
    Ok(result)
}

fn require_role_annotation(structure: &Structure) -> Result<(), HseError> {
    crate::chemistry::has_polymer_roles(structure)
        .then_some(())
        .ok_or(HseError::MissingRoleAnnotation)
}

fn role_atom(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
) -> Result<Option<AtomIndex>, HseError> {
    let mut matches = residue
        .atoms()
        .filter(|atom| {
            crate::chemistry::polymer_role(structure, atom.index().get())
                .is_some_and(|role| role.intersects(required))
        })
        .map(pdbiox_core::structure::AtomRef::index);
    let first = matches.next();
    if matches.next().is_some() {
        Err(HseError::AmbiguousRole {
            residue: residue.index(),
            role: required.code(),
        })
    } else {
        Ok(first)
    }
}

/// Adds neighbour `other` to the upper or lower count of residue `centre`.
fn classify(centres: &mut [Centre], centre: usize, other: usize) {
    let neighbour = centres[other].alpha;
    let host = &centres[centre];
    let to_neighbour = [
        neighbour[0] - host.alpha[0],
        neighbour[1] - host.alpha[1],
        neighbour[2] - host.alpha[2],
    ];
    let projection = to_neighbour[0] * host.direction[0]
        + to_neighbour[1] * host.direction[1]
        + to_neighbour[2] * host.direction[2];
    if projection > 0.0 {
        centres[centre].upper += 1;
    } else {
        centres[centre].lower += 1;
    }
}

#[cfg(test)]
#[path = "hse_tests.rs"]
mod tests;
