//! Radial distributions and coordination shells over explicit site selections.

use core::f64::consts::PI;
use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{PeriodicBox, SpatialBackend, SpatialError, pairs_within};

use crate::numeric::{f32_to_usize, f64_to_f32, u64_to_f64, usize_to_f32};

/// Explicit binning and normalization for a radial distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialDistributionOptions {
    /// Inclusive lower distance bound.
    pub minimum_distance: f32,
    /// Exclusive upper distance bound.
    pub maximum_distance: f32,
    /// Number of equal-width radial shells.
    pub bins: usize,
    /// Sample volume used to normalize shell counts into `g(r)`.
    pub volume: f64,
    /// Spatial implementation to use.
    pub backend: SpatialBackend,
}

/// Counts and normalized density for one spherical shell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialBin {
    /// Inclusive inner shell radius.
    pub lower: f32,
    /// Exclusive outer shell radius.
    pub upper: f32,
    /// Observed unique site pairs in this shell.
    pub count: u64,
    /// Pair density relative to an ideal uniform distribution.
    pub distribution: f64,
}

/// One explicitly ordered atom group represented by its centre of mass.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CentreGroup {
    /// Atoms contributing to the group centre.
    pub atoms: AtomSelection,
}

/// Why radial or coordination analysis could not be evaluated.
#[derive(Debug, thiserror::Error)]
pub enum RadialError {
    /// Distance bounds do not define a positive finite interval.
    #[error("radial bounds must be finite, non-negative, and increasing")]
    InvalidBounds,
    /// At least one shell is required.
    #[error("at least one radial bin is required")]
    NoBins,
    /// A finite positive normalization volume is required.
    #[error("normalization volume must be finite and positive")]
    InvalidVolume,
    /// The spatial query failed.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
    /// Groups or atom-aligned masses are empty, out of range or invalid.
    #[error("centre-of-mass groups and atom masses must be complete and physically valid")]
    InvalidGroups,
    /// The exact number of possible site pairs exceeds the supported integer domain.
    #[error("possible radial pair count exceeds the supported integer domain")]
    PairCountOverflow,
}

/// Computes a group-centre radial distribution using explicit atom masses.
///
/// Group coordinates must already be imaged as whole molecules. The resulting
/// centres can still be compared through the supplied periodic box.
///
/// # Errors
///
/// Returns invalid-group, radial-policy or spatial errors.
pub fn centre_of_mass_radial_distribution(
    positions: &[[f32; 3]],
    masses: &[f64],
    left: &[CentreGroup],
    right: &[CentreGroup],
    options: RadialDistributionOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<RadialBin>, RadialError> {
    if positions.len() != masses.len() || left.is_empty() || right.is_empty() {
        return Err(RadialError::InvalidGroups);
    }
    let same = left == right;
    let mut centres = group_centres(positions, masses, left)?;
    let left_selection =
        AtomSelection::All(u32::try_from(centres.len()).map_err(|_| RadialError::InvalidGroups)?);
    let right_selection = if same {
        left_selection.clone()
    } else {
        let start = u32::try_from(centres.len()).map_err(|_| RadialError::InvalidGroups)?;
        centres.extend(group_centres(positions, masses, right)?);
        let end = u32::try_from(centres.len()).map_err(|_| RadialError::InvalidGroups)?;
        AtomSelection::from_sorted((start..end).collect())
    };
    radial_distribution(
        &centres,
        &left_selection,
        &right_selection,
        options,
        periodic,
    )
}

fn group_centres(
    positions: &[[f32; 3]],
    masses: &[f64],
    groups: &[CentreGroup],
) -> Result<Vec<[f32; 3]>, RadialError> {
    groups
        .iter()
        .map(|group| {
            let atoms: Vec<_> = group.atoms.iter().collect();
            if atoms.is_empty() {
                return Err(RadialError::InvalidGroups);
            }
            let points = atoms
                .iter()
                .map(|atom| positions.get(*atom as usize).copied())
                .collect::<Option<Vec<_>>>()
                .ok_or(RadialError::InvalidGroups)?;
            let weights = atoms
                .iter()
                .map(|atom| masses.get(*atom as usize).copied())
                .collect::<Option<Vec<_>>>()
                .ok_or(RadialError::InvalidGroups)?;
            let centre =
                pdbiox_geom::centre_of_mass(&points, &weights).ok_or(RadialError::InvalidGroups)?;
            let centre = centre.map(f64_to_f32);
            centre
                .iter()
                .all(|value| value.is_finite())
                .then_some(centre)
                .ok_or(RadialError::InvalidGroups)
        })
        .collect()
}

/// Computes a normalized site-site radial distribution function.
///
/// Pair enumeration delegates to the shared spatial planner. When `left` and
/// `right` are identical, the ideal population is `N(N-1)/2`; otherwise it is
/// `N_left N_right`, matching the unique cross-pairs returned by the planner.
///
/// # Errors
///
/// Returns [`RadialError`] for invalid bounds, bin count, volume, selections,
/// or periodic/spatial input.
pub fn radial_distribution(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    options: RadialDistributionOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<RadialBin>, RadialError> {
    validate_options(options)?;
    let width = (options.maximum_distance - options.minimum_distance) / usize_to_f32(options.bins);
    let pairs = pairs_within(
        positions,
        left,
        right,
        options.maximum_distance,
        options.backend,
        periodic,
    )?;
    let mut counts = vec![0_u64; options.bins];

    for pair in pairs {
        let distance = pair.distance_squared.sqrt();
        if distance >= options.minimum_distance {
            let Some(bin) = f32_to_usize((distance - options.minimum_distance) / width) else {
                continue;
            };
            if let Some(count) = counts.get_mut(bin) {
                *count += 1;
            }
        }
    }

    let possible_pairs = u64_to_f64(possible_pair_count(left, right)?);
    Ok(counts
        .into_iter()
        .enumerate()
        .map(|(index, count)| make_bin(index, count, width, possible_pairs, options))
        .collect())
}

/// Counts neighbours in an explicit radial shell for every selected left site.
///
/// The output follows `left` selection order. For identical selections, each
/// unique pair contributes once to both endpoints.
///
/// # Errors
///
/// Returns [`RadialError`] for invalid shell bounds or spatial input.
pub fn coordination_numbers(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    minimum_distance: f32,
    maximum_distance: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<u32>, RadialError> {
    validate_bounds(minimum_distance, maximum_distance)?;
    let pairs = pairs_within(positions, left, right, maximum_distance, backend, periodic)?;
    let mut by_atom = vec![0_u32; positions.len()];

    for pair in pairs {
        if pair.distance_squared.sqrt() < minimum_distance {
            continue;
        }
        if left.contains(pair.first) && right.contains(pair.second) {
            by_atom[pair.first as usize] += 1;
        }
        if left.contains(pair.second) && right.contains(pair.first) {
            by_atom[pair.second as usize] += 1;
        }
    }

    Ok(left
        .into_iter()
        .map(|atom| by_atom[atom as usize])
        .collect())
}

fn validate_options(options: RadialDistributionOptions) -> Result<(), RadialError> {
    validate_bounds(options.minimum_distance, options.maximum_distance)?;
    if options.bins == 0 {
        return Err(RadialError::NoBins);
    }
    if !options.volume.is_finite() || options.volume <= 0.0 {
        return Err(RadialError::InvalidVolume);
    }
    Ok(())
}

fn validate_bounds(minimum: f32, maximum: f32) -> Result<(), RadialError> {
    if minimum.is_finite() && maximum.is_finite() && minimum >= 0.0 && maximum > minimum {
        Ok(())
    } else {
        Err(RadialError::InvalidBounds)
    }
}

fn possible_pair_count(left: &AtomSelection, right: &AtomSelection) -> Result<u64, RadialError> {
    let left_len = left.len();
    let right_len = right.len();
    if left == right {
        let predecessor = if left_len == 0 { 0 } else { left_len - 1 };
        left_len
            .checked_mul(predecessor)
            .map(|count| count / 2)
            .ok_or(RadialError::PairCountOverflow)
    } else {
        left_len
            .checked_mul(right_len)
            .ok_or(RadialError::PairCountOverflow)
    }
}

fn make_bin(
    index: usize,
    count: u64,
    width: f32,
    possible_pairs: f64,
    options: RadialDistributionOptions,
) -> RadialBin {
    let lower = options.minimum_distance + usize_to_f32(index) * width;
    let upper = lower + width;
    let shell_volume = 4.0 * PI * (f64::from(upper).powi(3) - f64::from(lower).powi(3)) / 3.0;
    let ideal = possible_pairs * shell_volume / options.volume;
    RadialBin {
        lower,
        upper,
        count,
        distribution: if ideal > 0.0 {
            u64_to_f64(count) / ideal
        } else {
            0.0
        },
    }
}

#[cfg(test)]
#[path = "radial_tests.rs"]
mod tests;
