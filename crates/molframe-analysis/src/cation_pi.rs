//! Cation–π interactions between annotated cations and aromatic rings.
//!
//! A positively charged group sitting over the face of an aromatic ring is
//! stabilised by the ring's π electrons. Cation membership and aromaticity come
//! from CCD-derived atom annotations and bonds. A pair counts only when
//! the cation is both close to the ring centre and roughly over its face —
//! judged by the angle between the ring normal and the line to the cation —
//! since a cation level with the ring edge is not a cation–π contact.

use molframe_core::index::ResidueIndex;
use molframe_core::structure::Structure;

use crate::numeric::f64_to_f32;
use std::collections::BTreeMap;

use crate::aromatic_stacking::{PiStackingError, aromatic_rings};

/// Explicit geometric policy for cation–π interactions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CationPiOptions {
    /// Maximum cation-to-ring-centre distance in ångström.
    pub maximum_distance: f32,
    /// Maximum angle from the ring normal in degrees.
    pub maximum_face_angle: f64,
    /// Numerical controls for aromatic ring plane fitting.
    pub plane_fit: molframe_geom::EigenOptions,
}

/// Why cation–π analysis could not be completed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum CationPiError {
    /// The requested geometry policy was invalid.
    #[error("cation-pi distance and face angle must be finite and geometrically valid")]
    InvalidOptions,
    /// Aromatic ring plane construction failed.
    #[error(transparent)]
    RingGeometry(#[from] PiStackingError),
}

/// A cation found over an aromatic ring.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CationPi {
    /// The residue carrying the cation.
    pub cation_residue: ResidueIndex,
    /// The aromatic residue.
    pub ring_residue: ResidueIndex,
    /// The distance from the cation to the ring centre, in ångström.
    pub distance: f32,
}

define_soa_table! {
    /// Native columnar storage for cation-pi interactions.
    pub struct CationPiTable for CationPi {
        /// Cation-bearing residue indices.
        cation_residue: ResidueIndex,
        /// Aromatic residue indices.
        ring_residue: ResidueIndex,
        /// Cation-to-ring distances.
        distance: f32,
    }
}

/// Finds cation–π interactions with the cation within `max_distance` of a ring.
///
/// Results are ordered by `(cation residue, ring residue)`. A cation off to the
/// side of a ring, rather than over its face, is not reported.
///
/// Runs in `O(cations · rings)` time, both of which are sparse.
///
/// # Errors
///
/// Returns [`CationPiError`] for invalid policy or failed ring-plane fitting.
pub fn cation_pi(
    structure: &Structure,
    options: CationPiOptions,
) -> Result<CationPiTable, CationPiError> {
    validate_options(options)?;
    let rings = aromatic_rings(structure, options.plane_fit)?;
    let limit = f64::from(options.maximum_distance);
    let mut found: BTreeMap<(ResidueIndex, ResidueIndex), f64> = BTreeMap::new();
    for atom in structure.data().atoms() {
        if crate::chemistry::formal_charge(structure, atom.index().get()).is_none_or(|v| v <= 0) {
            continue;
        }
        let (Some(residue), Some(position)) = (atom.residue(), atom.position()) else {
            continue;
        };
        let cation = position.map(f64::from);
        for ring in &rings {
            if ring.residue == residue.index() {
                continue;
            }
            let to_cation = [
                cation[0] - ring.centre[0],
                cation[1] - ring.centre[1],
                cation[2] - ring.centre[2],
            ];
            let distance = norm(to_cation);
            if distance == 0.0
                || distance > limit
                || off_axis_angle(to_cation, ring.normal, distance) > options.maximum_face_angle
            {
                continue;
            }
            found
                .entry((residue.index(), ring.residue))
                .and_modify(|current| *current = current.min(distance))
                .or_insert(distance);
        }
    }
    Ok(found
        .into_iter()
        .map(|((cation_residue, ring_residue), distance)| CationPi {
            cation_residue,
            ring_residue,
            distance: f64_to_f32(distance),
        })
        .collect())
}

fn validate_options(options: CationPiOptions) -> Result<(), CationPiError> {
    if options.maximum_distance.is_finite()
        && options.maximum_distance > 0.0
        && options.maximum_face_angle.is_finite()
        && (0.0..=90.0).contains(&options.maximum_face_angle)
    {
        Ok(())
    } else {
        Err(CationPiError::InvalidOptions)
    }
}

/// The angle in degrees between the ring axis and the line to the cation.
fn off_axis_angle(to_cation: [f64; 3], normal: [f64; 3], distance: f64) -> f64 {
    let dot = to_cation[0] * normal[0] + to_cation[1] * normal[1] + to_cation[2] * normal[2];
    let cosine = (dot / distance).abs().min(1.0);
    molframe_geom::degrees(cosine.acos())
}

/// Euclidean length of a vector.
fn norm(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

#[cfg(test)]
#[path = "cation_pi_tests.rs"]
mod tests;
