//! Isotropic B-factor consistency against explicit crystallographic TLS groups.

use super::BFactorError;
use pdbiox_core::index::AtomIndex;
use pdbiox_core::{AtomSelection, Structure};

/// Translation/libration/screw model in crystallographic TLS convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TlsModel {
    /// TLS origin in structure coordinate units.
    pub origin: [f64; 3],
    /// Symmetric translation tensor.
    pub translation: [[f64; 3]; 3],
    /// Symmetric libration tensor in radians squared.
    pub libration: [[f64; 3]; 3],
    /// Screw tensor.
    pub screw: [[f64; 3]; 3],
}

/// Explicit atoms governed by one TLS model.
#[derive(Clone, Debug, PartialEq)]
pub struct TlsGroup {
    /// Stable group identifier.
    pub id: String,
    /// Group atoms.
    pub atoms: AtomSelection,
    /// Declared TLS model.
    pub model: TlsModel,
}

/// One atom inconsistent with its declared TLS prediction.
#[derive(Clone, Debug, PartialEq)]
pub struct TlsBFactorFlag {
    /// TLS group identifier.
    pub group: String,
    /// Atom index.
    pub atom: AtomIndex,
    /// Recorded isotropic B factor.
    pub observed: f64,
    /// TLS-predicted isotropic B factor.
    pub predicted: f64,
    /// Signed observed-minus-predicted difference.
    pub deviation: f64,
}

/// TLS consistency findings and coverage.
#[derive(Clone, Debug, PartialEq)]
pub struct TlsBFactorReport {
    /// Group atom memberships intended for assessment.
    pub intended: usize,
    /// Memberships with coordinates and B factors.
    pub assessed: usize,
    /// Deviations beyond the explicit tolerance.
    pub flags: Vec<TlsBFactorFlag>,
}

/// Validates isotropic B factors against explicit TLS groups.
///
/// # Errors
///
/// Returns invalid matrices, coordinates or tolerance.
pub fn tls_b_factor_consistency(
    structure: &Structure,
    groups: &[TlsGroup],
    maximum_absolute_deviation: f64,
    symmetry_tolerance: f64,
) -> Result<TlsBFactorReport, BFactorError> {
    if !maximum_absolute_deviation.is_finite()
        || maximum_absolute_deviation < 0.0
        || !symmetry_tolerance.is_finite()
        || symmetry_tolerance < 0.0
        || groups
            .iter()
            .any(|group| !valid_model(group.model, symmetry_tolerance))
    {
        return Err(BFactorError::InvalidTls);
    }
    let mut intended = 0usize;
    let mut assessed = 0usize;
    let mut flags = Vec::new();
    for group in groups {
        for index in &group.atoms {
            intended += 1;
            let Some(atom) = structure.data().atom(AtomIndex::new(index)) else {
                continue;
            };
            let (Some(position), Some(observed)) = (atom.position(), atom.b_factor()) else {
                continue;
            };
            assessed += 1;
            let predicted = predicted_b(position.map(f64::from), group.model);
            let deviation = f64::from(observed) - predicted;
            if deviation.abs() > maximum_absolute_deviation {
                flags.push(TlsBFactorFlag {
                    group: group.id.clone(),
                    atom: atom.index(),
                    observed: f64::from(observed),
                    predicted,
                    deviation,
                });
            }
        }
    }
    Ok(TlsBFactorReport {
        intended,
        assessed,
        flags,
    })
}

fn predicted_b(position: [f64; 3], model: TlsModel) -> f64 {
    let relative: [f64; 3] = std::array::from_fn(|axis| position[axis] - model.origin[axis]);
    let [coordinate_x, coordinate_y, coordinate_z] = relative;
    let translation = model.translation;
    let libration = model.libration;
    let screw = model.screw;
    let trace = translation[0][0]
        + translation[1][1]
        + translation[2][2]
        + libration[0][0] * (coordinate_y * coordinate_y + coordinate_z * coordinate_z)
        + libration[1][1] * (coordinate_x * coordinate_x + coordinate_z * coordinate_z)
        + libration[2][2] * (coordinate_x * coordinate_x + coordinate_y * coordinate_y)
        - 2.0
            * (libration[0][1] * coordinate_x * coordinate_y
                + libration[0][2] * coordinate_x * coordinate_z
                + libration[1][2] * coordinate_y * coordinate_z)
        + 2.0
            * ((screw[1][0] - screw[0][1]) * coordinate_z
                + (screw[2][1] - screw[1][2]) * coordinate_x
                + (screw[0][2] - screw[2][0]) * coordinate_y);
    8.0 * core::f64::consts::PI.powi(2) * trace / 3.0
}

fn valid_model(model: TlsModel, tolerance: f64) -> bool {
    [
        model.origin.as_slice(),
        model.translation.as_flattened(),
        model.libration.as_flattened(),
        model.screw.as_flattened(),
    ]
    .into_iter()
    .flatten()
    .all(|value| value.is_finite())
        && symmetric(model.translation, tolerance)
        && symmetric(model.libration, tolerance)
}

fn symmetric(matrix: [[f64; 3]; 3], tolerance: f64) -> bool {
    (0..3).all(|row| {
        (row + 1..3).all(|column| (matrix[row][column] - matrix[column][row]).abs() <= tolerance)
    })
}
