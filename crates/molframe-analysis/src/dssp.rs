//! Secondary structure by the DSSP hydrogen-bond method.
//!
//! Backbone hydrogen bonds are found with the Kabsch–Sander electrostatic model:
//! a carbonyl C=O and an amide N–H bond when the energy of their four-atom
//! Coulomb interaction drops below −0.5 kcal/mol. The pattern of those bonds then
//! names the structure — a run of `i→i+4` bonds is an α-helix, reciprocal bonds
//! between distant residues are a β-bridge, and shorter turns are turns.
//!
//! The amide hydrogen is placed from the previous residue's carbonyl when the
//! model omits it. Bonds are searched within each chain, so the cost is
//! quadratic in a chain's length.

use molframe_chem::PolymerAtomRole;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::{ResidueRef, Structure};
use std::collections::BTreeSet;
use std::ops::RangeInclusive;

/// Explicit numerical and pattern definition of a DSSP-like assignment.
#[derive(Clone, Debug, PartialEq)]
pub struct DsspOptions {
    /// Electrostatic prefactor in energy-distance units.
    pub electrostatic_prefactor: f64,
    /// Interactions below this energy are hydrogen bonds.
    pub hydrogen_bond_energy: f64,
    /// Reconstructed amide N-H distance in angstroms.
    pub amide_hydrogen_distance: f32,
    /// Minimum residue-index separation for hydrogen-bond candidates.
    pub minimum_sequence_separation: usize,
    /// Donor/acceptor offset identifying the configured helix class.
    pub helix_offset: usize,
    /// Inclusive donor/acceptor offsets identifying turns.
    pub turn_offsets: RangeInclusive<usize>,
}

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

/// A residue's assigned secondary-structure state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SseKind {
    /// An α-helix residue.
    AlphaHelix,
    /// A β-strand (bridge) residue.
    Strand,
    /// A hydrogen-bonded turn.
    Turn,
    /// None of the above.
    Coil,
}

/// One residue and its secondary-structure state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SseRecord {
    /// The residue described.
    pub residue: ResidueIndex,
    /// Its assigned state.
    pub kind: SseKind,
}

define_soa_table! {
    /// Native columnar storage for secondary-structure assignments.
    pub struct SseTable for SseRecord {
        /// Residue indices.
        residue: ResidueIndex,
        /// Assigned secondary-structure states.
        kind: SseKind,
    }
}

/// Assigns secondary structure to every backbone residue of a structure.
///
/// Residues without the required semantic backbone roles are treated as coil.
/// Results are ordered by residue index.
///
/// Runs in `O(chain length²)` per chain.
///
/// # Errors
///
/// Returns [`DsspError::MissingRoleAnnotation`] when no explicit semantic role
/// profile has been applied, [`DsspError::InvalidOptions`] for an invalid
/// numerical definition, or [`DsspError::AmbiguousRole`] for non-unique
/// backbone roles.
pub fn secondary_structure(
    structure: &Structure,
    options: &DsspOptions,
) -> Result<SseTable, DsspError> {
    validate_options(options)?;
    if !crate::chemistry::has_polymer_roles(structure) {
        return Err(DsspError::MissingRoleAnnotation);
    }
    let mut records = Vec::new();
    for chain in structure.data().chains() {
        let residues: Vec<ResidueRef<'_>> = chain.residues().collect();
        let backbones = backbones(structure, &residues, options.amide_hydrogen_distance)?;
        let bonds = hydrogen_bonds(&backbones, options);
        let kinds = classify(&bonds, residues.len(), options);
        for (residue, kind) in residues.iter().zip(kinds) {
            records.push(SseRecord {
                residue: residue.index(),
                kind,
            });
        }
    }
    records.sort_by_key(|record| record.residue.get());
    Ok(records.into_iter().collect())
}

fn validate_options(options: &DsspOptions) -> Result<(), DsspError> {
    let valid = options.electrostatic_prefactor.is_finite()
        && options.electrostatic_prefactor > 0.0
        && options.hydrogen_bond_energy.is_finite()
        && options.amide_hydrogen_distance.is_finite()
        && options.amide_hydrogen_distance > 0.0
        && options.minimum_sequence_separation > 0
        && options.helix_offset > options.minimum_sequence_separation
        && !options.turn_offsets.is_empty()
        && *options.turn_offsets.start() > options.minimum_sequence_separation;
    valid.then_some(()).ok_or(DsspError::InvalidOptions)
}

/// The backbone atoms a residue needs, with the amide hydrogen estimated.
#[derive(Clone, Copy)]
struct Backbone {
    nitrogen: Option<[f32; 3]>,
    carbon: Option<[f32; 3]>,
    oxygen: Option<[f32; 3]>,
    hydrogen: Option<[f32; 3]>,
}

/// Extracts each residue's backbone, placing the amide H from the prior carbonyl.
fn backbones(
    structure: &Structure,
    residues: &[ResidueRef<'_>],
    hydrogen_distance: f32,
) -> Result<Vec<Backbone>, DsspError> {
    let mut backbones = Vec::with_capacity(residues.len());
    for (position, residue) in residues.iter().enumerate() {
        let nitrogen = role_position(structure, *residue, PolymerAtomRole::PROTEIN_NITROGEN)?;
        let previous_carbonyl = position
            .checked_sub(1)
            .and_then(|p| residues.get(p))
            .map(|previous| {
                Ok((
                    role_position(
                        structure,
                        *previous,
                        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                    )?,
                    role_position(
                        structure,
                        *previous,
                        PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                    )?,
                ))
            })
            .transpose()?
            .and_then(|(carbon, oxygen)| Some((carbon?, oxygen?)));
        let hydrogen = match (nitrogen, previous_carbonyl) {
            (Some(nitrogen), Some((carbon, oxygen))) => {
                Some(place_hydrogen(nitrogen, carbon, oxygen, hydrogen_distance))
            }
            _ => None,
        };
        backbones.push(Backbone {
            nitrogen,
            carbon: role_position(
                structure,
                *residue,
                PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
            )?,
            oxygen: role_position(
                structure,
                *residue,
                PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
            )?,
            hydrogen,
        });
    }
    Ok(backbones)
}

/// Places the amide hydrogen 1 Å from N, away from the previous carbonyl.
fn place_hydrogen(
    nitrogen: [f32; 3],
    carbon: [f32; 3],
    oxygen: [f32; 3],
    hydrogen_distance: f32,
) -> [f32; 3] {
    let direction = [
        carbon[0] - oxygen[0],
        carbon[1] - oxygen[1],
        carbon[2] - oxygen[2],
    ];
    let length =
        (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2])
            .sqrt();
    if length <= f32::EPSILON {
        return nitrogen;
    }
    [
        nitrogen[0] + hydrogen_distance * direction[0] / length,
        nitrogen[1] + hydrogen_distance * direction[1] / length,
        nitrogen[2] + hydrogen_distance * direction[2] / length,
    ]
}

/// The set of `(carbonyl residue, amide residue)` backbone hydrogen bonds.
fn hydrogen_bonds(backbones: &[Backbone], options: &DsspOptions) -> BTreeSet<(usize, usize)> {
    let mut bonds = BTreeSet::new();
    for donor in 0..backbones.len() {
        for acceptor in 0..backbones.len() {
            if donor.abs_diff(acceptor) < options.minimum_sequence_separation {
                continue;
            }
            if hbond_energy(
                &backbones[acceptor],
                &backbones[donor],
                options.electrostatic_prefactor,
            ) < options.hydrogen_bond_energy
            {
                bonds.insert((donor, acceptor));
            }
        }
    }
    bonds
}

/// The Kabsch–Sander energy between a carbonyl residue and an amide residue.
///
/// Returns a large positive value — no bond — when any of the four atoms is
/// missing.
fn hbond_energy(carbonyl: &Backbone, amide: &Backbone, prefactor: f64) -> f64 {
    let (Some(c), Some(o), Some(n), Some(h)) = (
        carbonyl.carbon,
        carbonyl.oxygen,
        amide.nitrogen,
        amide.hydrogen,
    ) else {
        return f64::INFINITY;
    };
    let on = molframe_geom::distance(o, n);
    let ch = molframe_geom::distance(c, h);
    let oh = molframe_geom::distance(o, h);
    let cn = molframe_geom::distance(c, n);
    if on <= 0.0 || ch <= 0.0 || oh <= 0.0 || cn <= 0.0 {
        return f64::INFINITY;
    }
    prefactor * (1.0 / on + 1.0 / ch - 1.0 / oh - 1.0 / cn)
}

/// Assigns a state to each residue from the hydrogen-bond pattern.
fn classify(bonds: &BTreeSet<(usize, usize)>, count: usize, options: &DsspOptions) -> Vec<SseKind> {
    let has = |i: usize, j: usize| bonds.contains(&(i, j));
    let mut kinds = vec![SseKind::Coil; count];

    // α-helix: residues bracketed by an i→i+4 backbone hydrogen bond.
    for &(i, j) in bonds {
        if j > i && j - i == options.helix_offset {
            for kind in &mut kinds[(i + 1)..j] {
                *kind = SseKind::AlphaHelix;
            }
        }
    }

    // β-bridges: reciprocal or offset bonds between residues over two apart.
    for i in 0..count {
        for j in 0..count {
            if i.abs_diff(j) <= 2 {
                continue;
            }
            let antiparallel = (has(i, j) && has(j, i))
                || (i >= 1 && j + 1 < count && has(i - 1, j + 1) && has(j - 1, i + 1));
            let parallel = (i >= 1 && has(i - 1, j) && has(j, i + 1))
                || (j >= 1 && has(j - 1, i) && has(i, j + 1));
            if antiparallel || parallel {
                mark_strand(&mut kinds, i);
                mark_strand(&mut kinds, j);
            }
        }
    }

    // Turns: residues bracketed by a shorter 3-, 4- or 5-bond, if still coil.
    for &(i, j) in bonds {
        if j > i && options.turn_offsets.contains(&(j - i)) {
            for kind in &mut kinds[(i + 1)..j] {
                if *kind == SseKind::Coil {
                    *kind = SseKind::Turn;
                }
            }
        }
    }
    kinds
}

/// Marks a residue as a strand unless it is already a helix.
fn mark_strand(kinds: &mut [SseKind], residue: usize) {
    if kinds[residue] != SseKind::AlphaHelix {
        kinds[residue] = SseKind::Strand;
    }
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
