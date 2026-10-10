//! Aromatic ring stacking between residues.
//!
//! Two aromatic rings interact when their centres come close and the rings take
//! one of two characteristic arrangements: stacked, with nearly parallel planes
//! lying one above the other (face-to-face or offset), or edge-to-face, the
//! "T-shaped" arrangement, with nearly perpendicular planes. Each ring is a
//! single cycle of the bond graph with its own centroid and plane normal, so a
//! tryptophan contributes two rings and a biaryl ligand two, and rings of one
//! residue are never paired with each other.
//!
//! Plane angle and centre distance alone do not identify a stack: two coplanar
//! rings side by side are parallel and close yet do not stack. The
//! classification therefore also looks at where the centroid-to-centroid vector
//! points relative to the planes:
//!
//! - **Parallel**: plane angle at most `maximum_parallel_angle`; the
//!   interplanar separation (the mean projection of the centroid vector on the
//!   two normals) lies within
//!   [`PiStackingGeometry::minimum_interplanar_separation`,
//!   `maximum_interplanar_separation`]; and the lateral offset (the centroid
//!   vector's component within the planes) is at most
//!   `maximum_lateral_offset`. Defaults are 3.0 to 4.5 ångström separation and
//!   3.5 ångström offset, which admit sandwich and displaced stacks and reject
//!   side-by-side coplanar rings.
//! - **T-shaped**: plane angle at least `minimum_perpendicular_angle`, and the
//!   centroid vector within `maximum_t_shaped_axis_angle` (default 40 degrees)
//!   of the normal of the ring whose face is pointed at, so the other ring's
//!   edge approaches that face. The pointed-at ring is the one whose normal
//!   lies closer to the centroid vector.
//!
//! Orientations in between are left unclassified rather than guessed. Rings are
//! compared directly; the cost is quadratic in the number of aromatic rings.
//!
//! Geometry is non-periodic by default. With
//! [`PiStackingGeometry::periodic`] the centroid vector uses the minimum image
//! of the structure's unit cell; rings are assumed to be whole in the
//! deposited coordinates.

use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::structure::Structure;
use molframe_spatial::PeriodicBox;

use crate::aromatic_rings::{AromaticRing, aromatic_rings};
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

/// Placement criteria that distinguish stacked and edge-to-face rings from
/// rings that merely have a compatible plane angle and a short distance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PiStackingGeometry {
    /// Smallest accepted distance between parallel ring planes, in ångström.
    /// Default 3.0.
    pub minimum_interplanar_separation: f64,
    /// Largest accepted distance between parallel ring planes, in ångström.
    /// Default 4.5.
    pub maximum_interplanar_separation: f64,
    /// Largest accepted in-plane displacement between parallel ring centres, in
    /// ångström. Default 3.5.
    pub maximum_lateral_offset: f64,
    /// Largest angle, in degrees, between the centroid vector and the normal of
    /// the ring presenting its face in a T-shaped pair. Default 40.
    pub maximum_t_shaped_axis_angle: f64,
    /// Apply the structure's unit cell and minimum-image convention. Default
    /// `false`; a missing or placeholder cell is then rejected as invalid options.
    pub periodic: bool,
}

impl Default for PiStackingGeometry {
    fn default() -> Self {
        Self {
            minimum_interplanar_separation: 3.0,
            maximum_interplanar_separation: 4.5,
            maximum_lateral_offset: 3.5,
            maximum_t_shaped_axis_angle: 40.0,
            periodic: false,
        }
    }
}

/// Why aromatic stacking analysis could not be completed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum PiStackingError {
    /// The requested geometry policy was invalid.
    #[error(
        "pi-stacking distances and angles must be finite, ordered and geometrically valid, \
         and periodic geometry needs a real unit cell"
    )]
    InvalidOptions,
    /// A ring plane could not be decomposed with the selected numerical profile.
    #[error("aromatic ring plane fitting failed: {0:?}")]
    Geometry(molframe_geom::EigenError),
}

/// The geometry of a ring-stacking interaction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StackingKind {
    /// The rings are stacked: nearly parallel planes, one above the other
    /// (face-to-face or offset).
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
    /// The lowest atom of the first ring, identifying it within its residue.
    pub first_atom: AtomIndex,
    /// The lowest atom of the second ring, identifying it within its residue.
    pub second_atom: AtomIndex,
    /// The distance between the ring centres, in ångström.
    pub centre_distance: f32,
    /// The angle between the ring planes, in degrees, from 0 to 90.
    pub angle: f64,
    /// Distance between the planes (parallel) or from the pointed-at ring's
    /// plane to the other centre (T-shaped), in ångström.
    pub interplanar_separation: f32,
    /// In-plane component of the centroid vector, in ångström.
    pub lateral_offset: f32,
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
        /// Lowest atom of each first ring.
        first_atom: AtomIndex,
        /// Lowest atom of each second ring.
        second_atom: AtomIndex,
        /// Ring-centre distances.
        centre_distance: f32,
        /// Acute plane angles.
        angle: f64,
        /// Interplanar separations.
        interplanar_separation: f32,
        /// Lateral offsets.
        lateral_offset: f32,
        /// Interaction classifications.
        kind: StackingKind,
    }
}

/// Finds aromatic ring stacks with centres within `maximum_centre_distance`,
/// classified with the default [`PiStackingGeometry`] (non-periodic).
///
/// Results are ordered by residue pair, then ring. Rings whose arrangement is
/// neither stacked nor edge-to-face within the recognised ranges are not
/// reported.
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
    pi_stacking_with_geometry(structure, options, PiStackingGeometry::default())
}

/// As [`pi_stacking`], with explicit placement criteria and periodicity.
///
/// # Errors
///
/// Returns [`PiStackingError`] for invalid policy, a missing or placeholder
/// cell when `geometry.periodic` is set, or failed ring-plane fitting.
pub fn pi_stacking_with_geometry(
    structure: &Structure,
    options: PiStackingOptions,
    geometry: PiStackingGeometry,
) -> Result<PiStackingTable, PiStackingError> {
    validate_options(options, geometry)?;
    let periodic = periodic_box(structure, geometry.periodic)?;
    let rings = aromatic_rings(structure, options.plane_fit).map_err(PiStackingError::Geometry)?;
    let limit = f64::from(options.maximum_centre_distance);
    let mut stacks = Vec::new();
    for (i, left) in rings.iter().enumerate() {
        for right in &rings[i + 1..] {
            if left.residue == right.residue {
                continue;
            }
            let (first, second) = if left.residue <= right.residue {
                (left, right)
            } else {
                (right, left)
            };
            let offset = displacement(periodic.as_ref(), first.centroid, second.centroid);
            let separation = norm(offset);
            if separation > limit || separation == 0.0 {
                continue;
            }
            let angle = plane_angle(first.normal, second.normal);
            let Some((kind, height)) =
                classify(first, second, offset, separation, angle, options, geometry)
            else {
                continue;
            };
            stacks.push(PiStacking {
                first: first.residue,
                second: second.residue,
                first_atom: first.first_atom(),
                second_atom: second.first_atom(),
                centre_distance: f64_to_f32(separation),
                angle,
                interplanar_separation: f64_to_f32(height),
                lateral_offset: f64_to_f32(
                    (separation * separation - height * height).max(0.0).sqrt(),
                ),
                kind,
            });
        }
    }
    stacks.sort_by_key(|stack| {
        (
            stack.first.get(),
            stack.second.get(),
            stack.first_atom.get(),
            stack.second_atom.get(),
        )
    });
    Ok(stacks.into_iter().collect())
}

/// Classifies one ring pair, returning the kind and its interplanar separation.
fn classify(
    first: &AromaticRing,
    second: &AromaticRing,
    offset: [f64; 3],
    separation: f64,
    angle: f64,
    options: PiStackingOptions,
    geometry: PiStackingGeometry,
) -> Option<(StackingKind, f64)> {
    let along_first = dot(offset, first.normal).abs();
    let along_second = dot(offset, second.normal).abs();
    if angle <= options.maximum_parallel_angle {
        let height = f64::midpoint(along_first, along_second);
        let lateral = (separation * separation - height * height).max(0.0).sqrt();
        let stacked = height >= geometry.minimum_interplanar_separation
            && height <= geometry.maximum_interplanar_separation
            && lateral <= geometry.maximum_lateral_offset;
        return stacked.then_some((StackingKind::Parallel, height));
    }
    if angle >= options.minimum_perpendicular_angle {
        let height = along_first.max(along_second);
        let axis_angle = molframe_geom::degrees((height / separation).min(1.0).acos());
        return (axis_angle <= geometry.maximum_t_shaped_axis_angle)
            .then_some((StackingKind::TShaped, height));
    }
    None
}

fn validate_options(
    options: PiStackingOptions,
    geometry: PiStackingGeometry,
) -> Result<(), PiStackingError> {
    if options.maximum_centre_distance.is_finite()
        && options.maximum_centre_distance > 0.0
        && options.maximum_parallel_angle.is_finite()
        && options.minimum_perpendicular_angle.is_finite()
        && options.maximum_parallel_angle >= 0.0
        && options.maximum_parallel_angle < options.minimum_perpendicular_angle
        && options.minimum_perpendicular_angle <= 90.0
        && geometry.minimum_interplanar_separation.is_finite()
        && geometry.minimum_interplanar_separation >= 0.0
        && geometry.maximum_interplanar_separation.is_finite()
        && geometry.minimum_interplanar_separation <= geometry.maximum_interplanar_separation
        && geometry.maximum_lateral_offset.is_finite()
        && geometry.maximum_lateral_offset >= 0.0
        && geometry.maximum_t_shaped_axis_angle.is_finite()
        && (0.0..=90.0).contains(&geometry.maximum_t_shaped_axis_angle)
    {
        Ok(())
    } else {
        Err(PiStackingError::InvalidOptions)
    }
}

/// The periodic box of a structure when periodic geometry is requested.
///
/// A missing cell and the placeholder unit cube some files write in place of a
/// cell are both rejected: neither describes real periodic images.
pub(crate) fn periodic_box(
    structure: &Structure,
    periodic: bool,
) -> Result<Option<PeriodicBox>, PiStackingError> {
    if !periodic {
        return Ok(None);
    }
    match structure.data().cell {
        Some(cell) if !cell.is_placeholder() => PeriodicBox::from_cell(cell)
            .map(Some)
            .map_err(|_| PiStackingError::InvalidOptions),
        _ => Err(PiStackingError::InvalidOptions),
    }
}

/// Vector from `from` to `to`, through the minimum image when periodic.
pub(crate) fn displacement(
    periodic: Option<&PeriodicBox>,
    from: [f64; 3],
    to: [f64; 3],
) -> [f64; 3] {
    match periodic {
        Some(cell) => cell.displacement_f64(from, to),
        None => [to[0] - from[0], to[1] - from[1], to[2] - from[2]],
    }
}

/// The acute angle in degrees between two plane normals.
fn plane_angle(first: [f64; 3], second: [f64; 3]) -> f64 {
    molframe_geom::degrees(dot(first, second).abs().min(1.0).acos())
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

#[cfg(test)]
#[path = "pi_stacking_tests.rs"]
mod tests;
