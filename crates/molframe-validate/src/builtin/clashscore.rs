//! Clashscore: serious steric overlaps per thousand atoms.
//!
//! A clash is a non-bonded pair whose van der Waals spheres overlap by at
//! least the tolerance (0.4 Å by default), found by the shared spatial clash
//! search, which already leaves out 1-2, 1-3 and 1-4 neighbours and atoms in
//! incompatible alternate locations.
//!
//! Two semantics, chosen by whether the model carries hydrogens:
//!
//! - **Heavy-atom**: no hydrogens are present (or they are ignored). Only
//!   heavy-atom pairs are counted and the score divides by the heavy-atom
//!   count. Hydrogen-bonded N/O pairs (including waters) sit at about 2.6-3.0
//!   Å, inside the tolerance of their summed radii, so a nitrogen or oxygen
//!   against another nitrogen or oxygen is exempt unless it overlaps by the
//!   tolerance plus `polar_margin`, which keeps a genuinely collapsed pair.
//! - **Hydrogen-aware**: hydrogens are present. Every pair is counted except a
//!   polar hydrogen against a nitrogen or oxygen, which is a hydrogen bond
//!   rather than a clash, and the score divides by all atoms.
//!
//! The two scores are not interchangeable and the result records which one
//! was computed, together with the options (the policy) that produced it.

use molframe_chem::RadiusSet;
use molframe_core::ExecutionContext;
use molframe_core::element::Element;
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;

use super::BuiltinError;
use crate::numeric::usize_to_f64;
use crate::steric_overlap::{ClashOptions, clashes_with};

/// Longest X-H bond used to find a hydrogen's parent atom, in ångström.
const PARENT_BOND: f64 = 1.3;

/// Whether hydrogens present in the model take part.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HydrogenMode {
    /// Use hydrogens when the model has any.
    Auto,
    /// Treat the model as heavy-atom only.
    Ignore,
}

/// Options for [`clash_score`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClashScoreOptions {
    /// Overlap at or above which a pair is a clash, in ångström.
    pub tolerance: f32,
    /// Radius set for the overlap.
    pub radius_set: RadiusSet,
    /// Hydrogen handling.
    pub hydrogens: HydrogenMode,
    /// Extra overlap beyond `tolerance` that a heavy-atom N/O pair needs to
    /// count, in ångström; unused when hydrogens take part.
    pub polar_margin: f32,
}

impl Default for ClashScoreOptions {
    fn default() -> Self {
        Self {
            tolerance: 0.4,
            radius_set: RadiusSet::Bondi,
            hydrogens: HydrogenMode::Auto,
            polar_margin: 0.6,
        }
    }
}

/// The clashscore and the residues behind it.
#[derive(Clone, Debug, PartialEq)]
pub struct ClashScore {
    /// Number of clashing pairs counted.
    pub clashes: usize,
    /// Atoms the score is normalised by.
    pub atoms: usize,
    /// Clashes per thousand atoms.
    pub score: f64,
    /// Whether hydrogens took part.
    pub hydrogen_aware: bool,
    /// Per residue, the number of counted clashes touching it.
    pub residue_clashes: Vec<u32>,
    /// Largest overlap among counted clashes, ångström.
    pub worst_overlap: f32,
    /// The policy the score was computed under.
    pub options: ClashScoreOptions,
}

struct Atoms {
    residue: Vec<u32>,
    element: Vec<Option<Element>>,
}

fn collect(structure: &Structure) -> Atoms {
    let count = structure.atom_count() as usize;
    let mut atoms = Atoms {
        residue: vec![u32::MAX; count],
        element: vec![None; count],
    };
    for residue in structure.data().residues() {
        for atom in residue.atoms() {
            let at = atom.index().as_usize();
            if let (Some(slot), Some(kind)) = (atoms.residue.get_mut(at), atoms.element.get_mut(at))
            {
                *slot = residue.index().get();
                *kind = atom.element();
            }
        }
    }
    atoms
}

fn is_polar(atoms: &Atoms, structure: &Structure, hydrogen: usize) -> bool {
    let positions = structure.positions();
    let Some(origin) = positions.get(hydrogen).copied() else {
        return false;
    };
    let Some(owner) = atoms.residue.get(hydrogen) else {
        return false;
    };
    let Some(residue) = structure
        .data()
        .residue(molframe_core::index::ResidueIndex::new(*owner))
    else {
        return false;
    };
    residue
        .atoms()
        .filter(|candidate| candidate.index().as_usize() != hydrogen)
        .filter_map(|candidate| {
            let position = positions.get(candidate.index().as_usize()).copied()?;
            let distance = molframe_geom::distance(origin, position);
            (distance <= PARENT_BOND).then_some((distance, candidate.element()))
        })
        .min_by(|left, right| left.0.total_cmp(&right.0))
        .is_some_and(|(_, element)| {
            element == Some(Element::NITROGEN) || element == Some(Element::OXYGEN)
        })
}

fn is_acceptor(kind: Option<Element>) -> bool {
    kind == Some(Element::NITROGEN) || kind == Some(Element::OXYGEN)
}

/// Computes the clashscore of a structure.
///
/// Runs in `O(atoms · local density)` time; the pair search is parallel
/// through the execution context.
///
/// # Errors
///
/// Returns [`BuiltinError`] if the neighbour search fails.
pub fn clash_score(
    structure: &Structure,
    options: ClashScoreOptions,
    context: &ExecutionContext,
) -> Result<ClashScore, BuiltinError> {
    let atoms = collect(structure);
    let hydrogen = |kind: &Option<Element>| kind.is_some_and(Element::is_hydrogen);
    let hydrogens = atoms.element.iter().filter(|kind| hydrogen(kind)).count();
    let aware = options.hydrogens == HydrogenMode::Auto && hydrogens > 0;
    let table = clashes_with(
        structure,
        ClashOptions::new(options.tolerance, options.radius_set, SpatialBackend::Auto),
        context,
    )?;
    let mut score = ClashScore {
        clashes: 0,
        score: 0.0,
        worst_overlap: 0.0,
        options,
        hydrogen_aware: aware,
        atoms: if aware {
            atoms.element.len()
        } else {
            atoms.element.len() - hydrogens
        },
        residue_clashes: vec![0; structure.residue_count()],
    };
    for clash in table.iter() {
        let (a, b) = (clash.first.as_usize(), clash.second.as_usize());
        let (kind_a, kind_b) = (
            atoms.element.get(a).copied().flatten(),
            atoms.element.get(b).copied().flatten(),
        );
        let (h_a, h_b) = (hydrogen(&kind_a), hydrogen(&kind_b));
        if (h_a || h_b) && !aware {
            continue;
        }
        let exempt = (h_a && is_acceptor(kind_b) && is_polar(&atoms, structure, a))
            || (h_b && is_acceptor(kind_a) && is_polar(&atoms, structure, b));
        let hydrogen_bond_pair = !aware
            && is_acceptor(kind_a)
            && is_acceptor(kind_b)
            && clash.overlap < options.tolerance + options.polar_margin;
        if exempt || hydrogen_bond_pair {
            continue;
        }
        score.clashes += 1;
        score.worst_overlap = score.worst_overlap.max(clash.overlap);
        let mut touched = [atoms.residue.get(a).copied(), atoms.residue.get(b).copied()];
        if touched[0] == touched[1] {
            touched[1] = None;
        }
        for residue in touched.into_iter().flatten() {
            if let Some(slot) = score.residue_clashes.get_mut(residue as usize) {
                *slot += 1;
            }
        }
    }
    if score.atoms > 0 {
        score.score = 1000.0 * usize_to_f64(score.clashes) / usize_to_f64(score.atoms);
    }
    Ok(score)
}

#[cfg(test)]
#[path = "clashscore_tests.rs"]
mod tests;
