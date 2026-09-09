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

use pdbiox_core::topology::{Csr, CsrBuilder};
use pdbiox_core::{ExecutionContext, selection::AtomSelection};
use pdbiox_spatial::{
    NeighborPair, SpatialBackend, SpatialError, SpatialSearchOptions,
    pairs_within_unsorted_with_options,
};

use crate::accessible_area::SasaError;
use crate::numeric::{f64_to_f32, usize_to_u32};

/// Probe-grown radii, atom centres, and the overlap adjacency they imply.
pub(crate) struct Neighbourhood<'a> {
    /// The stored coordinates, borrowed rather than promoted into a copy.
    ///
    /// A double-precision copy costs twenty-four bytes per atom on top of the
    /// twelve the coordinates already occupy, which at a billion atoms is more
    /// memory than the whole surface calculation needs. Promoting one centre
    /// when it is read is three conversions and touches no extra cache line
    /// (ADR-0025).
    positions: &'a [[f32; 3]],
    /// Each atom's radius grown by the probe, indexed by atom.
    pub(crate) expanded: Vec<f64>,
    /// Squared grown radii, kept beside `expanded` for hot point tests.
    pub(crate) expanded_squared: Vec<f64>,
    /// For each atom, the atoms whose grown spheres can reach it.
    pub(crate) adjacency: Csr<u32>,
}

impl Neighbourhood<'_> {
    /// One atom's centre in double precision.
    ///
    /// Out-of-range atoms return the origin, which no finite sphere contains
    /// at a positive radius, so a malformed index cannot report a false burial.
    pub(crate) fn centre(&self, atom: usize) -> [f64; 3] {
        match self.positions.get(atom) {
            Some(position) => promote(*position),
            None => [0.0; 3],
        }
    }

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
        for &neighbour in self.adjacency.row(atom) {
            let other = neighbour as usize;

            if excluded == Some(other) {
                continue;
            }

            if point_inside_sphere(point, self.centre(other), self.expanded_squared[other]) {
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
pub(crate) fn build<'a>(
    positions: &'a [[f32; 3]],
    radii: &[f32],
    probe: f32,
    context: &ExecutionContext,
) -> Result<Option<Neighbourhood<'a>>, SasaError> {
    validate_input(positions, radii, probe)?;

    if positions.is_empty() {
        return Ok(None);
    }

    let (expanded, widest) = expand_radii(radii, probe)?;
    let expanded_squared = expanded.iter().map(|radius| radius * radius).collect();

    let adjacency = neighbours(positions, &expanded, widest, context)?;

    Ok(Some(Neighbourhood {
        positions,
        expanded,
        expanded_squared,
        adjacency,
    }))
}

/// Validates neighbourhood input lengths and radii.
///
/// Runtime is `O(R)` and requires no allocation.
pub(crate) fn validate_input(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
) -> Result<(), SasaError> {
    if positions.len() != radii.len() {
        return Err(SasaError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
        });
    }

    validate_radii(radii, probe)
}

pub(crate) fn validate_radii(radii: &[f32], probe: f32) -> Result<(), SasaError> {
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
fn promote(position: [f32; 3]) -> [f64; 3] {
    [
        f64::from(position[0]),
        f64::from(position[1]),
        f64::from(position[2]),
    ]
}

/// Builds, for each atom, the atoms whose grown spheres can reach it.
///
/// A single search with the widest possible reach is a superset. Candidate
/// pairs are then filtered by their exact double-precision grown radii before
/// entering the adjacency, reducing all downstream local-density work.
fn neighbours(
    positions: &[[f32; 3]],
    expanded: &[f64],
    widest: f64,
    context: &ExecutionContext,
) -> Result<Csr<u32>, SpatialError> {
    if widest <= 0.0 {
        return Ok(CsrBuilder::with_degrees(&vec![0usize; positions.len()]).finish());
    }

    let all = AtomSelection::All(usize_to_u32(positions.len()));
    let cutoff = conservative_f32(2.0 * widest);

    // The pairs are scattered into per-atom adjacency and each list is sorted
    // afterwards, so the global pair order is never consumed: skip the sort.
    let pairs = pairs_within_unsorted_with_options(
        positions,
        &all,
        &all,
        cutoff,
        SpatialSearchOptions::with_backend(SpatialBackend::Auto),
        None,
        context,
    )?;

    // The counting pass sizes every row exactly, so the scatter below performs
    // no allocation at all: one flat buffer replaces a Vec header and a heap
    // allocation per atom.
    let degrees = overlap_degrees(&pairs, positions, expanded);
    let mut builder = CsrBuilder::with_degrees(&degrees);

    fill_adjacency(&mut builder, pairs, positions, expanded);

    Ok(builder.finish_canonical())
}

/// Counts exact sphere-overlap degrees to reserve adjacency capacity.
///
/// Runtime is `O(P)` for spatial candidate pairs and output space is `O(A)`.
fn overlap_degrees(pairs: &[NeighborPair], positions: &[[f32; 3]], expanded: &[f64]) -> Vec<usize> {
    let mut degrees = vec![0usize; positions.len()];

    for pair in pairs {
        if !pair_overlaps(pair, positions, expanded) {
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
/// row capacities.
fn fill_adjacency(
    adjacency: &mut CsrBuilder<u32>,
    pairs: Vec<NeighborPair>,
    positions: &[[f32; 3]],
    expanded: &[f64],
) {
    for pair in pairs {
        if !pair_overlaps(&pair, positions, expanded) {
            continue;
        }

        adjacency.push(pair.first as usize, pair.second);
        adjacency.push(pair.second as usize, pair.first);
    }
}

/// Returns whether the grown spheres associated with `pair` genuinely overlap.
///
/// The exact test is performed in double precision. Tangential spheres need not
/// enter the adjacency because they occlude no finite surface patch.
fn pair_overlaps(pair: &NeighborPair, positions: &[[f32; 3]], expanded: &[f64]) -> bool {
    let first = pair.first as usize;
    let second = pair.second as usize;

    let (Some(&first_point), Some(&second_point), Some(&first_radius), Some(&second_radius)) = (
        positions.get(first),
        positions.get(second),
        expanded.get(first),
        expanded.get(second),
    ) else {
        return false;
    };

    let reach = first_radius + second_radius;

    squared_distance(promote(first_point), promote(second_point)) < reach * reach
}

/// Tests whether `point` lies strictly inside a sphere.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(crate) fn point_inside_sphere(point: [f64; 3], centre: [f64; 3], radius_squared: f64) -> bool {
    squared_distance(point, centre) < radius_squared
}

/// Computes squared Euclidean distance between double-precision points.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(crate) fn squared_distance(left: [f64; 3], right: [f64; 3]) -> f64 {
    let dx = left[0] - right[0];
    let dy = left[1] - right[1];
    let dz = left[2] - right[2];

    dx * dx + dy * dy + dz * dz
}

/// Converts a non-negative `f64` bound to a conservative `f32` upper bound.
///
/// If normal rounding moves downward, the next representable positive `f32` is
/// selected so the broad-phase neighbour search cannot lose a boundary pair.
pub(crate) fn conservative_f32(value: f64) -> f32 {
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
