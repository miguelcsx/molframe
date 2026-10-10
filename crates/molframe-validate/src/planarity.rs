//! Aromatic ring planarity.
//!
//! CCD-declared aromatic systems are flat by chemistry; a ring that has
//! puckered in a model is almost always a refinement artefact. Each aromatic
//! ring — a cycle of the bond graph, not the whole aromatic part of a residue —
//! is fitted to its own best plane and flagged when the root-mean-square
//! departure from that plane exceeds a tolerance. Rings of one residue are
//! therefore judged independently: a biaryl twisted about its linking bond has
//! two flat rings and is not flagged, while a tryptophan whose pyrrole ring has
//! buckled is flagged without being masked by the benzene ring.
//!
//! A structure without bonds has no ring topology and yields no flags.

use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::structure::Structure;

/// An aromatic ring that has departed from planarity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarityFlag {
    /// The residue owning the non-planar ring.
    pub residue: ResidueIndex,
    /// The lowest atom of the non-planar ring, identifying it within the residue.
    pub first_atom: AtomIndex,
    /// The RMS distance of the ring atoms from their best plane, in ångström.
    pub deviation: f64,
}

/// Explicit policy for aromatic planarity validation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarityOptions {
    /// Maximum accepted RMS departure from the fitted plane, in ångström.
    pub maximum_deviation: f64,
    /// Numerical controls for plane fitting.
    pub plane_fit: molframe_geom::EigenOptions,
}

/// Why aromatic planarity validation could not be completed.
#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
pub enum PlanarityError {
    /// The distance threshold was not a finite, non-negative value.
    #[error("planarity tolerance must be finite and non-negative")]
    InvalidTolerance,
    /// A ring plane could not be decomposed.
    #[error("aromatic ring plane fitting failed: {0:?}")]
    Geometry(molframe_geom::EigenError),
}

/// Flags aromatic rings that are non-planar beyond `maximum_deviation`.
///
/// A tolerance of a few hundredths of an ångström catches real puckering while
/// tolerating coordinate noise. Results are ordered by residue index, then ring.
///
/// Runs in time linear in the number of aromatic bonds.
///
/// # Errors
///
/// Returns [`PlanarityError`] for an invalid tolerance or failed plane fit.
pub fn nonplanar_aromatic_rings(
    structure: &Structure,
    options: PlanarityOptions,
) -> Result<Vec<PlanarityFlag>, PlanarityError> {
    if !options.maximum_deviation.is_finite() || options.maximum_deviation < 0.0 {
        return Err(PlanarityError::InvalidTolerance);
    }
    let mut flags = Vec::new();
    let rings = molframe_analysis::aromatic_rings(structure, options.plane_fit)
        .map_err(PlanarityError::Geometry)?;
    for ring in rings {
        let points: Option<Vec<[f32; 3]>> = ring
            .atoms
            .iter()
            .map(|atom| {
                structure
                    .data()
                    .atom(*atom)
                    .and_then(molframe_core::structure::AtomRef::position)
            })
            .collect();
        let Some(points) = points else {
            continue;
        };
        let Some(deviation) =
            molframe_geom::plane_deviation_with_options(&points, options.plane_fit)
                .map_err(PlanarityError::Geometry)?
        else {
            continue;
        };
        if deviation > options.maximum_deviation {
            flags.push(PlanarityFlag {
                residue: ring.residue,
                first_atom: ring.first_atom(),
                deviation,
            });
        }
    }
    flags.sort_by_key(|flag| (flag.residue.get(), flag.first_atom.get()));
    Ok(flags)
}

#[cfg(test)]
#[path = "planarity_tests.rs"]
mod tests;
