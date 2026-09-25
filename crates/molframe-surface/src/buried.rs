//! Surface buried when two groups of atoms come together.
//!
//! The interface area is what the solvent stops seeing once a complex forms:
//! the accessible area of each partner measured alone, less the accessible area
//! of the two measured together. Because the two partners each occlude part of
//! the other, that difference is positive wherever they touch and zero when they
//! do not.
//!
//! By convention the buried surface is reported as the total area removed from
//! both partners, not the area of one face, so it is the full sum of the two
//! contributions rather than half of it.
//!
//! Cost is three solvent-accessible surface evaluations, so it inherits the
//! `O(atoms · points · local density)` cost of the underlying construction.

use crate::accessible_area::{SasaError, shrake_rupley};
use molframe_core::ExecutionContext;

/// Explicit participation of an atom in a two-molecule buried-surface calculation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoleculeRole {
    /// Atom belongs to the first molecule.
    First,
    /// Atom belongs to the second molecule.
    Second,
    /// Atom is deliberately excluded from all surface evaluations.
    Excluded,
}

/// Why a role-partitioned buried surface could not be computed.
#[derive(Debug, thiserror::Error)]
pub enum BuriedSurfaceError {
    /// Positions, radii, and molecule roles do not describe the same atoms.
    #[error(
        "expected one radius and role for each of {positions} positions, found {radii} radii and {roles} roles"
    )]
    LengthMismatch {
        /// Number of supplied positions.
        positions: usize,
        /// Number of supplied radii.
        radii: usize,
        /// Number of supplied roles.
        roles: usize,
    },
    /// One of the underlying surface evaluations failed.
    #[error(transparent)]
    Surface(#[from] SasaError),
}

/// The accessible areas that decide an interface, and the area it buries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BuriedSurface {
    /// Accessible area of the first group measured in isolation.
    pub first_alone: f64,
    /// Accessible area of the second group measured in isolation.
    pub second_alone: f64,
    /// Accessible area of both groups measured together.
    pub together: f64,
    /// Area removed from the solvent by contact: `first_alone + second_alone - together`.
    pub buried: f64,
}

/// Computes the surface buried between the atoms marked `true` and the rest.
///
/// `in_first` selects the first group position by position; every atom not
/// marked joins the second group. An atom in neither a crowded interface nor a
/// clash contributes nothing, so a partition with an empty side buries no area.
///
/// Runs in three accessible-surface passes over the atoms.
///
/// # Errors
///
/// Returns [`SasaError::LengthMismatch`] when the partition, positions and radii
/// disagree on length, and otherwise the same errors as [`shrake_rupley`].
pub fn buried_surface(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    points: u16,
    in_first: &[bool],
    context: &ExecutionContext,
) -> Result<BuriedSurface, SasaError> {
    if positions.len() != radii.len() || positions.len() != in_first.len() {
        return Err(SasaError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
        });
    }

    let mut first_positions = Vec::new();
    let mut first_radii = Vec::new();
    let mut second_positions = Vec::new();
    let mut second_radii = Vec::new();
    for (index, &first) in in_first.iter().enumerate() {
        if first {
            first_positions.push(positions[index]);
            first_radii.push(radii[index]);
        } else {
            second_positions.push(positions[index]);
            second_radii.push(radii[index]);
        }
    }

    let first_alone = total(&first_positions, &first_radii, probe, points, context)?;
    let second_alone = total(&second_positions, &second_radii, probe, points, context)?;
    let together = total(positions, radii, probe, points, context)?;
    let buried = first_alone + second_alone - together;

    Ok(BuriedSurface {
        first_alone,
        second_alone,
        together,
        buried,
    })
}

/// Computes solvent-excluded surface buried between two explicitly labelled molecules.
///
/// The partners are evaluated independently and as a complex. Atoms labelled
/// [`MoleculeRole::Excluded`] take part in none of the evaluations, making
/// waters, ions, membranes, and environment atoms an explicit caller choice.
///
/// # Errors
///
/// Returns [`BuriedSurfaceError::LengthMismatch`] when atom-aligned inputs differ
/// in length, and otherwise forwards errors from
/// [`solvent_excluded_surface`](crate::solvent_excluded_surface).
pub fn buried_solvent_excluded_surface(
    positions: &[[f32; 3]],
    radii: &[f32],
    roles: &[MoleculeRole],
    probe: f32,
    resolution: f32,
) -> Result<BuriedSurface, BuriedSurfaceError> {
    buried_solvent_excluded_surface_with_options(
        positions,
        radii,
        roles,
        probe,
        crate::SurfaceGridOptions::standard(resolution),
    )
}

/// Computes buried solvent-excluded area with an explicit grid resource budget.
///
/// # Errors
///
/// Returns [`BuriedSurfaceError::LengthMismatch`] for atom-aligned input drift
/// and otherwise forwards strict surface-grid errors.
pub fn buried_solvent_excluded_surface_with_options(
    positions: &[[f32; 3]],
    radii: &[f32],
    roles: &[MoleculeRole],
    probe: f32,
    grid_options: crate::SurfaceGridOptions,
) -> Result<BuriedSurface, BuriedSurfaceError> {
    if positions.len() != radii.len() || positions.len() != roles.len() {
        return Err(BuriedSurfaceError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
            roles: roles.len(),
        });
    }

    let (first_positions, first_radii) = select_role(positions, radii, roles, MoleculeRole::First);
    let (second_positions, second_radii) =
        select_role(positions, radii, roles, MoleculeRole::Second);
    let (together_positions, together_radii) = select_participants(positions, radii, roles);

    let first_alone = crate::solvent_excluded_surface_with_options(
        &first_positions,
        &first_radii,
        probe,
        grid_options,
    )?
    .area;
    let second_alone = crate::solvent_excluded_surface_with_options(
        &second_positions,
        &second_radii,
        probe,
        grid_options,
    )?
    .area;
    let together = crate::solvent_excluded_surface_with_options(
        &together_positions,
        &together_radii,
        probe,
        grid_options,
    )?
    .area;

    Ok(BuriedSurface {
        first_alone,
        second_alone,
        together,
        buried: first_alone + second_alone - together,
    })
}

fn select_role(
    positions: &[[f32; 3]],
    radii: &[f32],
    roles: &[MoleculeRole],
    selected: MoleculeRole,
) -> (Vec<[f32; 3]>, Vec<f32>) {
    positions
        .iter()
        .copied()
        .zip(radii.iter().copied())
        .zip(roles.iter().copied())
        .filter_map(|((position, radius), role)| (role == selected).then_some((position, radius)))
        .unzip()
}

fn select_participants(
    positions: &[[f32; 3]],
    radii: &[f32],
    roles: &[MoleculeRole],
) -> (Vec<[f32; 3]>, Vec<f32>) {
    positions
        .iter()
        .copied()
        .zip(radii.iter().copied())
        .zip(roles.iter().copied())
        .filter_map(|((position, radius), role)| {
            (role != MoleculeRole::Excluded).then_some((position, radius))
        })
        .unzip()
}

/// Sums the per-atom accessible areas of one group.
fn total(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    points: u16,
    context: &ExecutionContext,
) -> Result<f64, SasaError> {
    let areas = shrake_rupley(positions, radii, probe, points, context)?;
    Ok(areas.iter().sum())
}

#[cfg(test)]
#[path = "buried_tests.rs"]
mod tests;
