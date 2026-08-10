//! Planarity validation for caller- or CCD-supplied explicit atom planes.

use crate::{PlanarityError, PlanarityOptions};
use pdbiox_core::{AtomSelection, Structure};

/// One explicit plane restraint, independent of aromaticity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlaneRestraint {
    /// Stable caller identifier.
    pub id: String,
    /// Atoms intended to share a plane.
    pub atoms: AtomSelection,
}

/// One assessed plane outside its explicit tolerance.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneRestraintFlag {
    /// Stable caller identifier.
    pub id: String,
    /// RMS distance from the fitted plane.
    pub deviation: f64,
}

/// Findings and restraint coverage.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneRestraintReport {
    /// Supplied restraints.
    pub intended: usize,
    /// Restraints with at least three positioned atoms.
    pub assessed: usize,
    /// Outliers in supplied order.
    pub flags: Vec<PlaneRestraintFlag>,
}

/// Validates any explicit plane restraints without aromatic-name inference.
///
/// # Errors
///
/// Returns invalid tolerance or eigensolver errors.
pub fn plane_restraint_outliers(
    structure: &Structure,
    restraints: &[PlaneRestraint],
    options: PlanarityOptions,
) -> Result<PlaneRestraintReport, PlanarityError> {
    if !options.maximum_deviation.is_finite() || options.maximum_deviation < 0.0 {
        return Err(PlanarityError::InvalidTolerance);
    }
    let mut assessed = 0usize;
    let mut flags = Vec::new();
    for restraint in restraints {
        let points: Vec<_> = restraint
            .atoms
            .iter()
            .filter_map(|index| {
                structure
                    .data()
                    .atom(pdbiox_core::index::AtomIndex::new(index))?
                    .position()
            })
            .collect();
        let Some(deviation) = pdbiox_geom::plane_deviation_with_options(&points, options.plane_fit)
            .map_err(PlanarityError::Geometry)?
        else {
            continue;
        };
        assessed += 1;
        if deviation > options.maximum_deviation {
            flags.push(PlaneRestraintFlag {
                id: restraint.id.clone(),
                deviation,
            });
        }
    }
    Ok(PlaneRestraintReport {
        intended: restraints.len(),
        assessed,
        flags,
    })
}
