//! Aromatic ring stacking between residues.
//!
//! Two aromatic rings interact when their centres come close and their planes
//! take one of two characteristic orientations: nearly parallel, for the
//! face-to-face and offset stacks, or nearly perpendicular, for the edge-to-face
//! "T-shaped" arrangement. Each ring is reduced to a centre and a plane normal,
//! and a pair is reported only when both the distance and the angle fall in a
//! recognised range — orientations in between are left unclassified rather than
//! guessed.
//!
//! Aromatic residues are sparse, so the ring pairs are compared directly; the
//! cost is quadratic in the number of aromatic rings, not in the atom count.

use molframe_core::index::ResidueIndex;
use molframe_core::structure::Structure;

use crate::numeric::f64_to_f32;

/// Explicit geometric policy for aromatic ring stacking.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PiStackingOptions {
    /// Maximum ring-centre separation in ångström.
    pub maximum_centre_distance: f32,
    /// Largest acute plane angle classified as parallel.
    pub maximum_parallel_angle: f64,
    /// Smallest acute plane angle classified as T-shaped.
    pub minimum_perpendicular_angle: f64,
    /// Numerical controls for aromatic ring plane fitting.
    pub plane_fit: molframe_geom::EigenOptions,
}

/// Why aromatic stacking analysis could not be completed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum PiStackingError {
    /// The requested geometry policy was invalid.
    #[error("pi-stacking distances and angles must be finite, ordered and geometrically valid")]
    InvalidOptions,
    /// A ring plane could not be decomposed with the selected numerical profile.
    #[error("aromatic ring plane fitting failed: {0:?}")]
    Geometry(molframe_geom::EigenError),
}

/// The geometry of a ring-stacking interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StackingKind {
    /// The ring planes are nearly parallel (face-to-face or offset).
    Parallel,
    /// The ring planes are nearly perpendicular (edge-to-face).
    TShaped,
}

/// Two aromatic rings found stacking.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PiStacking {
    /// The lower-indexed residue.
    pub first: ResidueIndex,
    /// The higher-indexed residue.
    pub second: ResidueIndex,
    /// The distance between the ring centres, in ångström.
    pub centre_distance: f32,
    /// The angle between the ring planes, in degrees, from 0 to 90.
    pub angle: f64,
    /// Which arrangement it is.
    pub kind: StackingKind,
}

define_soa_table! {
    /// Native columnar storage for aromatic stacking interactions.
    pub struct PiStackingTable for PiStacking {
        /// Lower residue indices.
        first: ResidueIndex,
        /// Higher residue indices.
        second: ResidueIndex,
        /// Ring-centre distances.
        centre_distance: f32,
        /// Acute plane angles.
        angle: f64,
        /// Interaction classifications.
        kind: StackingKind,
    }
}

/// One aromatic ring reduced to its centre and plane normal.
pub(crate) struct Ring {
    pub(crate) residue: ResidueIndex,
    pub(crate) centre: [f64; 3],
    pub(crate) normal: [f64; 3],
}

/// Finds aromatic ring stacks with centres within `max_centre_distance`.
///
/// Results are ordered by residue pair. Rings whose planes are neither parallel
/// nor perpendicular within the recognised ranges are not reported.
///
/// Runs in `O(rings²)` time in the number of aromatic rings.
///
/// # Errors
///
/// Returns [`PiStackingError`] for invalid policy or failed ring-plane fitting.
pub fn pi_stacking(
    structure: &Structure,
    options: PiStackingOptions,
) -> Result<PiStackingTable, PiStackingError> {
    validate_options(options)?;
    let rings = aromatic_rings(structure, options.plane_fit)?;
    let limit = f64::from(options.maximum_centre_distance);
    let mut stacks = Vec::new();
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            let separation = distance(rings[i].centre, rings[j].centre);
            if separation > limit {
                continue;
            }
            let angle = plane_angle(rings[i].normal, rings[j].normal);
            let kind = if angle <= options.maximum_parallel_angle {
                StackingKind::Parallel
            } else if angle >= options.minimum_perpendicular_angle {
                StackingKind::TShaped
            } else {
                continue;
            };
            stacks.push(PiStacking {
                first: rings[i].residue,
                second: rings[j].residue,
                centre_distance: f64_to_f32(separation),
                angle,
                kind,
            });
        }
    }
    stacks.sort_by_key(|stack| (stack.first.get(), stack.second.get()));
    Ok(stacks.into_iter().collect())
}

fn validate_options(options: PiStackingOptions) -> Result<(), PiStackingError> {
    if options.maximum_centre_distance.is_finite()
        && options.maximum_centre_distance > 0.0
        && options.maximum_parallel_angle.is_finite()
        && options.minimum_perpendicular_angle.is_finite()
        && options.maximum_parallel_angle >= 0.0
        && options.maximum_parallel_angle < options.minimum_perpendicular_angle
        && options.minimum_perpendicular_angle <= 90.0
    {
        Ok(())
    } else {
        Err(PiStackingError::InvalidOptions)
    }
}

/// Reduces every aromatic residue with a fittable ring to a centre and normal.
pub(crate) fn aromatic_rings(
    structure: &Structure,
    plane_fit: molframe_geom::EigenOptions,
) -> Result<Vec<Ring>, PiStackingError> {
    let mut rings = Vec::new();
    for residue in structure.data().residues() {
        let points: Vec<_> = residue
            .atoms()
            .filter(|atom| crate::chemistry::is_aromatic(structure, atom.index().get()))
            .filter_map(molframe_core::structure::AtomRef::position)
            .collect();
        let Some(plane) = molframe_geom::best_fit_plane_with_options(&points, plane_fit)
            .map_err(PiStackingError::Geometry)?
        else {
            continue;
        };
        rings.push(Ring {
            residue: residue.index(),
            centre: plane.centre,
            normal: plane.normal,
        });
    }
    Ok(rings)
}

/// The acute angle in degrees between two plane normals.
fn plane_angle(first: [f64; 3], second: [f64; 3]) -> f64 {
    let dot = first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
    molframe_geom::degrees(dot.abs().min(1.0).acos())
}

/// Euclidean distance between two points.
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
#[path = "pi_stacking_tests.rs"]
mod tests;
