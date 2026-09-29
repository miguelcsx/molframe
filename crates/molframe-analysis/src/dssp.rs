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
use num_traits::ToPrimitive;
use std::collections::BTreeSet;
use std::ops::RangeInclusive;

const CA_GRID_CUTOFF: f32 = 9.0;

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
    /// The residue could not be evaluated from the available backbone data.
    Unknown,
    /// An α-helix residue.
    AlphaHelix,
    /// A β-strand (bridge) residue.
    Strand,
    /// A hydrogen-bonded turn.
    Turn,
    /// None of the above.
    Coil,
}

impl SseKind {
    /// Whether this value represents an evaluated assignment rather than a missing backbone result.
    #[must_use]
    pub const fn is_evaluated(self) -> bool {
        !matches!(self, Self::Unknown)
    }
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

impl SseTable {
    /// Returns the computed DSSP rows as typed placement input.
    /// This reuses the existing result and performs no coordinate or frame work.
    pub fn placements(&self) -> impl Iterator<Item = (ResidueIndex, SseKind)> + '_ {
        self.residue()
            .iter()
            .copied()
            .zip(self.kind().iter().copied())
    }
}

/// Assigns DSSP-like secondary structure to every residue in `structure`.
///
/// The assignment uses the Kabsch–Sander electrostatic energy of backbone
/// carbonyl and amide groups, followed by the hydrogen-bond pattern rules in
/// `options`. Backbone roles are resolved from the structure's explicit
/// semantic polymer-role annotations; an amide hydrogen is estimated from the
/// preceding carbonyl when the structure does not provide one.
///
/// Residues are processed independently within each chain and the returned
/// rows preserve chain/residue traversal order. A residue whose required
/// backbone roles are not all present is reported as [`SseKind::Unknown`]. An
/// evaluable residue that matches no helix, strand, or turn is reported as
/// [`SseKind::Coil`].
///
/// The hydrogen-bond search takes `O(n²)` time for a chain containing `n`
/// residues. The caller must provide options with finite, meaningful numeric
/// values; use [`DsspOptions`] to start from the standard definition.
///
/// # Errors
///
/// Returns [`DsspError::InvalidOptions`] if any numerical option is non-finite
/// or outside its allowed domain. Returns
/// [`DsspError::MissingRoleAnnotation`] when the structure has no explicit
/// semantic polymer atom-role profile. Returns [`DsspError::AmbiguousRole`]
/// when a residue contains more than one atom for a backbone role that must be
/// unique.
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
        let evaluable = backbones
            .iter()
            .map(Backbone::is_evaluable)
            .collect::<Vec<_>>();
        let kinds = classify(&bonds, &evaluable, options);
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
    ca: Option<[f32; 3]>,
    nitrogen: Option<[f32; 3]>,
    carbon: Option<[f32; 3]>,
    oxygen: Option<[f32; 3]>,
    hydrogen: Option<[f32; 3]>,
}

impl Backbone {
    const fn is_evaluable(&self) -> bool {
        self.nitrogen.is_some()
            && self.carbon.is_some()
            && self.oxygen.is_some()
            && self.hydrogen.is_some()
    }
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
            ca: role_position(structure, *residue, PolymerAtomRole::PROTEIN_ALPHA_CARBON)?,
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
    let mut entries = backbones
        .iter()
        .enumerate()
        .filter_map(|(residue, backbone)| {
            let position = backbone.ca?;
            let cell = ca_cell(position)?;
            Some(CaEntry { cell, residue })
        })
        .collect::<Vec<_>>();
    entries.sort_unstable();

    let mut bonds = BTreeSet::new();
    for left in &entries {
        let Some(left_backbone) = backbones.get(left.residue) else {
            continue;
        };
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let cell = [
                        left.cell[0].saturating_add(dx),
                        left.cell[1].saturating_add(dy),
                        left.cell[2].saturating_add(dz),
                    ];
                    let range = ca_range(&entries, cell);
                    for right in &entries[range] {
                        if right.residue <= left.residue
                            || right.residue.abs_diff(left.residue)
                                < options.minimum_sequence_separation
                        {
                            continue;
                        }
                        let Some(right_backbone) = backbones.get(right.residue) else {
                            continue;
                        };
                        if ca_distance_sq(left_backbone.ca, right_backbone.ca)
                            > CA_GRID_CUTOFF * CA_GRID_CUTOFF
                        {
                            continue;
                        }
                        if hbond_energy(
                            left_backbone,
                            right_backbone,
                            options.electrostatic_prefactor,
                        ) < options.hydrogen_bond_energy
                        {
                            bonds.insert((right.residue, left.residue));
                        }
                        if hbond_energy(
                            right_backbone,
                            left_backbone,
                            options.electrostatic_prefactor,
                        ) < options.hydrogen_bond_energy
                        {
                            bonds.insert((left.residue, right.residue));
                        }
                    }
                }
            }
        }
    }
    bonds
}

/// The Kabsch–Sander energy between a carbonyl residue and an amide residue.
///
/// Returns a large positive value — no bond — when any of the four atoms is
/// missing.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct CaEntry {
    cell: [i32; 3],
    residue: usize,
}

fn ca_range(entries: &[CaEntry], cell: [i32; 3]) -> std::ops::Range<usize> {
    let start = entries.partition_point(|entry| entry.cell < cell);
    let end = entries.partition_point(|entry| entry.cell <= cell);
    start..end
}

fn ca_cell(position: [f32; 3]) -> Option<[i32; 3]> {
    if !position.iter().all(|value| value.is_finite()) {
        return None;
    }
    Some([
        (position[0] / CA_GRID_CUTOFF).floor().to_i32()?,
        (position[1] / CA_GRID_CUTOFF).floor().to_i32()?,
        (position[2] / CA_GRID_CUTOFF).floor().to_i32()?,
    ])
}

fn ca_distance_sq(left: Option<[f32; 3]>, right: Option<[f32; 3]>) -> f32 {
    let (Some(left), Some(right)) = (left, right) else {
        return f32::INFINITY;
    };
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}

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
fn classify(
    bonds: &BTreeSet<(usize, usize)>,
    evaluable: &[bool],
    options: &DsspOptions,
) -> Vec<SseKind> {
    let count = evaluable.len();
    let has = |i: usize, j: usize| bonds.contains(&(i, j));
    let mut kinds = evaluable
        .iter()
        .map(|value| {
            if *value {
                SseKind::Coil
            } else {
                SseKind::Unknown
            }
        })
        .collect::<Vec<_>>();

    // α-helix: residues bracketed by an i→i+4 backbone hydrogen bond.
    for &(i, j) in bonds {
        if j > i && j - i == options.helix_offset {
            for kind in &mut kinds[(i + 1)..j] {
                *kind = SseKind::AlphaHelix;
            }
        }
    }

    // β-bridges: reciprocal or offset bonds between residues over two apart.
    // Iterate the sparse bond set rather than the full residue Cartesian product.
    for &(i, j) in bonds {
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
    if kinds[residue] != SseKind::AlphaHelix && kinds[residue] != SseKind::Unknown {
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
