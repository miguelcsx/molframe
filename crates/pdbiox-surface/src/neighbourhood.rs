//! The geometry every surface computation shares: probe-grown radii, atom
//! centres in double precision, and the per-atom list of atoms whose grown
//! spheres can reach them.
//!
//! Growing each radius by the probe once, and finding overlaps once, is work
//! that the accessible-surface sampler, the Lee–Richards integrator and the
//! molecular-surface pass all need identically. Keeping it here means a single
//! definition of "which atoms can shade this one" rather than three that could
//! drift apart.
//!
//! The neighbour list comes from the shared spatial search rather than an
//! all-pairs scan, so building it costs `O(atoms · local density)` rather than
//! quadratic in the atom count.

use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{NeighborPair, SpatialBackend, SpatialError, pairs_within};

use crate::accessible_area::SasaError;
use crate::numeric::{f64_to_f32, usize_to_u32};

/// Probe-grown radii, atom centres, and the overlap adjacency they imply.
pub(crate) struct Neighbourhood {
    /// Atom centres promoted to double precision, indexed by atom.
    pub(crate) centres: Vec<[f64; 3]>,
    /// Each atom's radius grown by the probe, indexed by atom.
    pub(crate) expanded: Vec<f64>,
    /// For each atom, the atoms whose grown spheres can reach it.
    pub(crate) adjacency: Vec<Vec<u32>>,
}

impl Neighbourhood {
    /// Whether `point` lies outside every one of `atom`'s neighbours' grown
    /// spheres, and so is not buried by any of them.
    ///
    /// Distances are compared squared against the squared grown radius, so no
    /// square root is taken per neighbour.
    pub(crate) fn point_is_clear(&self, point: [f64; 3], atom: usize) -> bool {
        self.point_is_clear_except(point, atom, None)
    }

    /// Performs the same test while deliberately omitting one neighbour.
    ///
    /// Runtime is `O(D)` for the atom's overlap degree and performs no
    /// allocation.
    pub(crate) fn point_is_clear_except(
        &self,
        point: [f64; 3],
        atom: usize,
        excluded: Option<usize>,
    ) -> bool {
        for &neighbour in &self.adjacency[atom] {
            let other = neighbour as usize;

            if excluded == Some(other) {
                continue;
            }

            if point_inside_sphere(point, self.centres[other], self.expanded[other]) {
                return false;
            }
        }

        true
    }
}

/// Validates the input and builds the shared neighbourhood.
///
/// Returns `Ok(None)` for an empty atom set — valid, but with nothing to
/// measure. The radii are the bare atomic radii; the probe is added here, so a
/// zero probe leaves the van der Waals spheres untouched.
///
/// # Errors
///
/// Returns [`SasaError::LengthMismatch`] when the counts disagree,
/// [`SasaError::InvalidProbe`] for a non-finite or negative probe,
/// [`SasaError::InvalidRadius`] for a non-finite or negative radius, and
/// [`SasaError::Spatial`] when the neighbour search rejects the input.
pub(crate) fn build(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
) -> Result<Option<Neighbourhood>, SasaError> {
    validate_input(positions, radii, probe)?;

    if positions.is_empty() {
        return Ok(None);
    }

    let (expanded, widest) = expand_radii(radii, probe)?;
    let centres = promote_centres(positions);

    let adjacency = neighbours(positions, &centres, &expanded, widest)?;

    Ok(Some(Neighbourhood {
        centres,
        expanded,
        adjacency,
    }))
}

/// Validates neighbourhood input lengths and radii.
///
/// Runtime is `O(R)` and requires no allocation.
fn validate_input(positions: &[[f32; 3]], radii: &[f32], probe: f32) -> Result<(), SasaError> {
    if positions.len() != radii.len() {
        return Err(SasaError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
        });
    }

    if !probe.is_finite() || probe < 0.0 {
        return Err(SasaError::InvalidProbe);
    }

    if radii
        .iter()
        .any(|radius| !radius.is_finite() || *radius < 0.0)
    {
        return Err(SasaError::InvalidRadius);
    }

    Ok(())
}

/// Grows all atomic radii and returns the largest grown radius.
///
/// Runtime and output space are `O(R)`.
fn expand_radii(radii: &[f32], probe: f32) -> Result<(Vec<f64>, f64), SasaError> {
    let probe = f64::from(probe);
    let mut expanded = Vec::with_capacity(radii.len());
    let mut widest = 0.0f64;

    for &radius in radii {
        if !radius.is_finite() || radius < 0.0 {
            return Err(SasaError::InvalidRadius);
        }

        let grown = f64::from(radius) + probe;
        widest = widest.max(grown);
        expanded.push(grown);
    }

    Ok((expanded, widest))
}

/// Promotes atom centres to double precision once for all surface algorithms.
///
/// Runtime and output space are `O(A)`.
fn promote_centres(positions: &[[f32; 3]]) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|position| {
            [
                f64::from(position[0]),
                f64::from(position[1]),
                f64::from(position[2]),
            ]
        })
        .collect()
}

/// Builds, for each atom, the atoms whose grown spheres can reach it.
///
/// A single search with the widest possible reach is a superset. Candidate
/// pairs are then filtered by their exact double-precision grown radii before
/// entering the adjacency, reducing all downstream local-density work.
fn neighbours(
    positions: &[[f32; 3]],
    centres: &[[f64; 3]],
    expanded: &[f64],
    widest: f64,
) -> Result<Vec<Vec<u32>>, SpatialError> {
    let mut adjacency = vec![Vec::new(); positions.len()];

    if widest <= 0.0 {
        return Ok(adjacency);
    }

    let all = AtomSelection::All(usize_to_u32(positions.len()));
    let cutoff = conservative_f32(2.0 * widest);

    let pairs = pairs_within(positions, &all, &all, cutoff, SpatialBackend::Auto, None)?;

    let degrees = overlap_degrees(&pairs, centres, expanded);

    adjacency = degrees.into_iter().map(Vec::with_capacity).collect();

    fill_adjacency(&mut adjacency, pairs, centres, expanded);

    canonicalise_adjacency(&mut adjacency);

    Ok(adjacency)
}

/// Counts exact sphere-overlap degrees to reserve adjacency capacity.
///
/// Runtime is `O(P)` for spatial candidate pairs and output space is `O(A)`.
fn overlap_degrees(pairs: &[NeighborPair], centres: &[[f64; 3]], expanded: &[f64]) -> Vec<usize> {
    let mut degrees = vec![0usize; centres.len()];

    for pair in pairs {
        if !pair_overlaps(pair, centres, expanded) {
            continue;
        }

        if let Some(degree) = degrees.get_mut(pair.first as usize) {
            *degree += 1;
        }

        if let Some(degree) = degrees.get_mut(pair.second as usize) {
            *degree += 1;
        }
    }

    degrees
}

/// Populates exact overlap adjacency from spatial candidate pairs.
///
/// Runtime is `O(P)` and no allocation occurs beyond the previously reserved
/// adjacency capacities.
fn fill_adjacency(
    adjacency: &mut [Vec<u32>],
    pairs: Vec<NeighborPair>,
    centres: &[[f64; 3]],
    expanded: &[f64],
) {
    for pair in pairs {
        if !pair_overlaps(&pair, centres, expanded) {
            continue;
        }

        if let Some(list) = adjacency.get_mut(pair.first as usize) {
            list.push(pair.second);
        }

        if let Some(list) = adjacency.get_mut(pair.second as usize) {
            list.push(pair.first);
        }
    }
}

/// Sorts and deduplicates every per-atom neighbour list.
///
/// This defensively removes repeated pair orientations without changing
/// accessibility semantics.
fn canonicalise_adjacency(adjacency: &mut [Vec<u32>]) {
    for neighbours in adjacency {
        if neighbours.len() < 2 {
            continue;
        }

        neighbours.sort_unstable();
        neighbours.dedup();
    }
}

/// Returns whether the grown spheres associated with `pair` genuinely overlap.
///
/// The exact test is performed in double precision. Tangential spheres need not
/// enter the adjacency because they occlude no finite surface patch.
fn pair_overlaps(pair: &NeighborPair, centres: &[[f64; 3]], expanded: &[f64]) -> bool {
    let first = pair.first as usize;
    let second = pair.second as usize;

    let (Some(&first_centre), Some(&second_centre), Some(&first_radius), Some(&second_radius)) = (
        centres.get(first),
        centres.get(second),
        expanded.get(first),
        expanded.get(second),
    ) else {
        return false;
    };

    let reach = first_radius + second_radius;

    squared_distance(first_centre, second_centre) < reach * reach
}

/// Tests whether `point` lies strictly inside a sphere.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn point_inside_sphere(point: [f64; 3], centre: [f64; 3], radius: f64) -> bool {
    squared_distance(point, centre) < radius * radius
}

/// Computes squared Euclidean distance between double-precision points.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn squared_distance(left: [f64; 3], right: [f64; 3]) -> f64 {
    let dx = left[0] - right[0];
    let dy = left[1] - right[1];
    let dz = left[2] - right[2];

    dx * dx + dy * dy + dz * dz
}

/// Converts a non-negative `f64` bound to a conservative `f32` upper bound.
///
/// If normal rounding moves downward, the next representable positive `f32` is
/// selected so the broad-phase neighbour search cannot lose a boundary pair.
fn conservative_f32(value: f64) -> f32 {
    let rounded = f64_to_f32(value);

    if !rounded.is_finite() || f64::from(rounded) >= value || rounded <= 0.0 {
        return rounded;
    }

    if rounded.to_bits() == f32::MAX.to_bits() {
        f32::INFINITY
    } else {
        f32::from_bits(rounded.to_bits() + 1)
    }
}
