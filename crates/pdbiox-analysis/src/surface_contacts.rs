//! Contacts supported by exposed surface near the interatomic axis.

use crate::Contact;
use pdbiox_core::structure::Structure;
use pdbiox_core::{ExecutionContext, selection::AtomSelection};
use pdbiox_spatial::{
    PairQuery, SpatialBackend, SpatialSearchOptions, reduce_pairs_within_unsorted,
};
use pdbiox_surface::{ExcludedSurfacePoint, SasaError, surface_points_excluding_pairs};
use std::collections::BTreeMap;

use crate::numeric::{f32_to_usize, f64_to_f32, usize_to_u32};

/// Controls exposed-surface contact detection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceContactOptions {
    /// Extra separation accepted beyond the two atomic radii.
    pub tolerance: f32,
    /// Solvent probe radius.
    pub probe: f32,
    /// Surface samples per square ångström.
    pub surface_density: f32,
    /// Minimum exposed patch area on both atoms.
    pub minimum_area: f32,
    /// Spatial query implementation.
    pub backend: SpatialBackend,
}

/// Finds contacts that satisfy both a radius cutoff and exposed-surface test.
///
/// `min_area` is the minimum exposed patch on each atom and `surface_density`
/// is measured in points per square ångström. A candidate survives only when
/// both atoms have the rounded number of points within `sqrt(min_area / π)` of
/// the segment joining their centres.
///
/// # Errors
///
/// Returns [`SasaError`] for invalid radii, probe, density or spatial inputs.
pub fn surface_contacts(
    structure: &Structure,
    radii: &[f32],
    options: SurfaceContactOptions,
    context: &ExecutionContext,
) -> Result<Vec<Contact>, SasaError> {
    if radii.len() != structure.atom_count() as usize {
        return Err(SasaError::LengthMismatch {
            positions: structure.atom_count() as usize,
            radii: radii.len(),
        });
    }
    if !options.tolerance.is_finite()
        || options.tolerance < 0.0
        || !options.minimum_area.is_finite()
        || options.minimum_area <= 0.0
    {
        return Err(SasaError::InvalidRadius);
    }
    let positions = structure.positions();
    let maximum_radius = radii.iter().copied().fold(0.0_f32, f32::max);
    let maximum_cutoff = 2.0 * maximum_radius + options.tolerance;
    let all = AtomSelection::All(structure.atom_count());

    // Only pairs whose spheres actually touch survive, a small fraction of the
    // candidates. Per-block vectors are merged in stable block order.
    let parts = reduce_pairs_within_unsorted(
        &PairQuery {
            positions,
            left: &all,
            right: &all,
            cutoff: maximum_cutoff,
            options: SpatialSearchOptions::with_backend(options.backend),
            periodic: None,
            context,
        },
        Vec::new,
        |close, pair| {
            let first = pair.first as usize;
            let second = pair.second as usize;
            let (Some(&a), Some(&b)) = (positions.get(first), positions.get(second)) else {
                return;
            };
            let distance = f64_to_f32(pdbiox_geom::distance(a, b));
            if distance <= radii[first] + radii[second] + options.tolerance {
                close.push((first, second, distance));
            }
        },
    )?;
    let mut close: Vec<_> = parts.into_iter().flatten().collect();

    // The streaming query leaves pairs unordered; ordering the survivors keeps
    // the excluded-point sampling and the reported contacts deterministic.
    close.sort_by_key(|&(first, second, _)| (first, second));
    let pair_indices: Vec<_> = close
        .iter()
        .map(|&(first, second, _)| (first, second))
        .collect();
    let points = surface_points_excluding_pairs(
        positions,
        radii,
        options.probe,
        options.surface_density,
        &pair_indices,
        context,
    )?;
    let mut by_pair: BTreeMap<(usize, usize), Vec<ExcludedSurfacePoint>> = BTreeMap::new();
    for point in points {
        by_pair
            .entry((point.atom, point.excluded))
            .or_default()
            .push(point);
    }
    let max_delta_squared = options.minimum_area / core::f32::consts::PI;
    let required = f32_to_usize(
        (options.surface_density * options.minimum_area)
            .round()
            .max(1.0),
    )
    .ok_or(SasaError::InvalidDensity)?;
    let mut contacts = Vec::new();
    for (first, second, distance) in close {
        let a = positions[first];
        let b = positions[second];
        if support(
            by_pair.get(&(first, second)).map(Vec::as_slice),
            a,
            b,
            max_delta_squared,
            required,
        ) && support(
            by_pair.get(&(second, first)).map(Vec::as_slice),
            a,
            b,
            max_delta_squared,
            required,
        ) {
            contacts.push(Contact {
                first: pdbiox_core::index::AtomIndex::new(usize_to_u32(first)),
                second: pdbiox_core::index::AtomIndex::new(usize_to_u32(second)),
                distance,
            });
        }
    }
    Ok(contacts)
}

fn support(
    points: Option<&[ExcludedSurfacePoint]>,
    start: [f32; 3],
    end: [f32; 3],
    max_delta_squared: f32,
    required: usize,
) -> bool {
    let Some(points) = points else {
        return false;
    };
    let mut count = 0usize;
    for point in points {
        if point_segment_distance_squared(point.position, start, end) <= max_delta_squared {
            count += 1;
            if count >= required {
                return true;
            }
        }
    }
    false
}

fn point_segment_distance_squared(point: [f32; 3], start: [f32; 3], end: [f32; 3]) -> f32 {
    let axis = [end[0] - start[0], end[1] - start[1], end[2] - start[2]];
    let relative = [
        point[0] - start[0],
        point[1] - start[1],
        point[2] - start[2],
    ];
    let length_squared = axis.iter().map(|value| value * value).sum::<f32>();
    let parameter = if length_squared <= f32::EPSILON {
        0.0
    } else {
        relative
            .iter()
            .zip(axis)
            .map(|(left, right)| left * right)
            .sum::<f32>()
            / length_squared
    }
    .clamp(0.0, 1.0);
    (0..3)
        .map(|dimension| {
            let delta = point[dimension] - (start[dimension] + parameter * axis[dimension]);
            delta * delta
        })
        .sum()
}

#[cfg(test)]
#[path = "surface_contacts_tests.rs"]
mod tests;
