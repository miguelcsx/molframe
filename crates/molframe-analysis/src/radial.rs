//! Radial distributions and coordination shells over explicit site selections.

use core::f64::consts::PI;
use molframe_core::{ExecutionContext, selection::AtomSelection};
use molframe_spatial::{
    PairQuery, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions,
    for_each_pairs_within_unsorted,
};
use std::collections::BTreeMap;

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

/// Distance shell and backend for coordination counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoordinationOptions {
    /// Inclusive lower shell bound.
    pub minimum_distance: f32,
    /// Inclusive upper shell bound.
    pub maximum_distance: f32,
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
/// Groups are identified by their atom set. A group that appears in both lists,
/// or more than once in one list, is a single site: it is never paired with
/// itself, a pair of sites is counted once however many lists contain them, and
/// the ideal pair population follows the same unique-pair rule as
/// [`radial_distribution`].
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
    context: &ExecutionContext,
) -> Result<Vec<RadialBin>, RadialError> {
    if positions.len() != masses.len() || left.is_empty() || right.is_empty() {
        return Err(RadialError::InvalidGroups);
    }
    let mut sites: Vec<&CentreGroup> = Vec::new();
    let mut identity: BTreeMap<Vec<u32>, u32> = BTreeMap::new();
    let left_ids = site_ids(left, &mut sites, &mut identity)?;
    let right_ids = site_ids(right, &mut sites, &mut identity)?;
    let mut centres = Vec::with_capacity(sites.len());
    for group in sites {
        centres.extend(group_centres(
            positions,
            masses,
            core::slice::from_ref(group),
        )?);
    }
    radial_distribution(
        &centres,
        &AtomSelection::from_sorted(left_ids),
        &AtomSelection::from_sorted(right_ids),
        options,
        periodic,
        context,
    )
}

/// Assigns each group a site id by atom-set identity and returns the sorted,
/// unique ids of the given list.
fn site_ids<'a>(
    groups: &'a [CentreGroup],
    sites: &mut Vec<&'a CentreGroup>,
    identity: &mut BTreeMap<Vec<u32>, u32>,
) -> Result<Vec<u32>, RadialError> {
    let mut ids = Vec::with_capacity(groups.len());
    for group in groups {
        let key: Vec<u32> = group.atoms.iter().collect();
        let next = u32::try_from(sites.len()).map_err(|_| RadialError::InvalidGroups)?;
        let id = *identity.entry(key).or_insert(next);
        if id == next {
            sites.push(group);
        }
        ids.push(id);
    }
    ids.sort_unstable();
    ids.dedup();
    Ok(ids)
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
            let centre = molframe_geom::centre_of_mass(&points, &weights)
                .ok_or(RadialError::InvalidGroups)?;
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
/// `right` are identical, the ideal population is `N(N-1)/2`; for disjoint
/// selections it is `N_left N_right`. Overlapping selections share `I` sites:
/// the planner never pairs a site with itself and reports a pair of shared
/// sites once, so the ideal population is `N_left N_right - I - I(I-1)/2`.
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
    context: &ExecutionContext,
) -> Result<Vec<RadialBin>, RadialError> {
    validate_options(options)?;
    let mut counts = vec![0_u64; options.bins];
    let width = (options.maximum_distance - options.minimum_distance) / usize_to_f32(options.bins);
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions,
            left,
            right,
            cutoff: options.maximum_distance,
            options: SpatialSearchOptions::with_backend(options.backend),
            periodic,
            context,
        },
        |pair| {
            let distance = pair.distance_squared.sqrt();
            if distance >= options.minimum_distance {
                let Some(bin) = f32_to_usize((distance - options.minimum_distance) / width) else {
                    return;
                };
                if let Some(count) = counts.get_mut(bin) {
                    *count += 1;
                }
            }
        },
    )?;

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
    options: CoordinationOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Vec<u32>, RadialError> {
    validate_bounds(options.minimum_distance, options.maximum_distance)?;
    let minimum_squared = options.minimum_distance * options.minimum_distance;
    let same_selection = left == right;
    let mut by_atom = vec![0_u32; positions.len()];
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions,
            left,
            right,
            cutoff: options.maximum_distance,
            options: SpatialSearchOptions::with_backend(options.backend),
            periodic,
            context,
        },
        |pair| {
            if pair.distance_squared < minimum_squared {
                return;
            }
            if same_selection {
                by_atom[pair.first as usize] += 1;
                by_atom[pair.second as usize] += 1;
                return;
            }
            if left.contains(pair.first) && right.contains(pair.second) {
                by_atom[pair.first as usize] += 1;
            }
            if left.contains(pair.second) && right.contains(pair.first) {
                by_atom[pair.second as usize] += 1;
            }
        },
    )?;

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

/// Counts the unique unordered pairs of distinct sites with one end in `left`
/// and the other in `right`.
///
/// Ordered pairs `(l, r)` with `l != r` number `|L||R| - |I|` for an
/// intersection `I`. A pair of two distinct shared sites appears as both
/// `(a, b)` and `(b, a)`, so `C(|I|, 2)` of them are removed.
fn possible_pair_count(left: &AtomSelection, right: &AtomSelection) -> Result<u64, RadialError> {
    let shared = left.intersect(right).len();
    let ordered = left
        .len()
        .checked_mul(right.len())
        .ok_or(RadialError::PairCountOverflow)?;
    let shared_pairs = shared
        .checked_mul(shared.saturating_sub(1))
        .ok_or(RadialError::PairCountOverflow)?
        / 2;
    Ok(ordered - shared - shared_pairs)
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
