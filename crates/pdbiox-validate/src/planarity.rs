//! Aromatic ring planarity.
//!
//! CCD-declared aromatic systems are flat by chemistry; a ring that has
//! puckered in a model is almost always a refinement
//! artefact. Each ring's atoms are fitted to their best plane and the residue is
//! flagged when the root-mean-square departure from that plane exceeds a
//! tolerance.
//!
//! A ring with fewer than three of its atoms present cannot be fitted and is
//! skipped rather than guessed at.

use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::Structure;
use pdbiox_core::{AtomAnnotation, Presence};

/// An aromatic residue whose ring has departed from planarity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarityFlag {
    /// The residue whose ring is non-planar.
    pub residue: ResidueIndex,
    /// The RMS distance of the ring atoms from their best plane, in ångström.
    pub deviation: f64,
}

/// Explicit policy for aromatic planarity validation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarityOptions {
    /// Maximum accepted RMS departure from the fitted plane, in ångström.
    pub maximum_deviation: f64,
    /// Numerical controls for plane fitting.
    pub plane_fit: pdbiox_geom::EigenOptions,
}

/// Why aromatic planarity validation could not be completed.
#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
pub enum PlanarityError {
    /// The distance threshold was not a finite, non-negative value.
    #[error("planarity tolerance must be finite and non-negative")]
    InvalidTolerance,
    /// A ring plane could not be decomposed.
    #[error("aromatic ring plane fitting failed: {0:?}")]
    Geometry(pdbiox_geom::EigenError),
}

/// Flags aromatic residues whose ring is non-planar beyond `tolerance`.
///
/// A tolerance of a few hundredths of an ångström catches real puckering while
/// tolerating coordinate noise. Results are ordered by residue index.
///
/// Runs in `O(residues)` time.
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
    let aromatic = structure
        .annotations()
        .get(pdbiox_core::AROMATIC_ATOM_ANNOTATION);
    for residue in structure.data().residues() {
        let points: Vec<_> = residue
            .atoms()
            .filter(|atom| aromatic_atom(aromatic, atom.index().get()))
            .filter_map(pdbiox_core::structure::AtomRef::position)
            .collect();
        let Some(deviation) = pdbiox_geom::plane_deviation_with_options(&points, options.plane_fit)
            .map_err(PlanarityError::Geometry)?
        else {
            continue;
        };
        if deviation > options.maximum_deviation {
            flags.push(PlanarityFlag {
                residue: residue.index(),
                deviation,
            });
        }
    }
    flags.sort_by_key(|flag| flag.residue.get());
    Ok(flags)
}

fn aromatic_atom(annotation: Option<&AtomAnnotation>, atom: u32) -> bool {
    let Some(AtomAnnotation::Boolean(column)) = annotation else {
        return false;
    };
    matches!(column.get(atom), Some((true, Presence::Present)))
}

#[cfg(test)]
#[path = "planarity_tests.rs"]
mod tests;
