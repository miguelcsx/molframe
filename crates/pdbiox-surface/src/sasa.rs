//! Solvent-accessible surface area by the Shrake–Rupley construction.
//!
//! Each atom is grown by the probe radius and dusted with a fixed set of test
//! points on that expanded sphere. A point counts as accessible when it lies
//! outside every neighbour's expanded sphere; the accessible fraction times the
//! sphere area is the atom's contribution. Summing the per-atom areas gives the
//! surface of the whole set.
//!
//! Neighbours come from the shared spatial search rather than an all-pairs scan,
//! so the cost is `O(atoms · points · local density)` rather than quadratic in
//! the atom count. The double-precision point tests keep the result stable
//! regardless of the order neighbours arrive in.

use pdbiox_spatial::SpatialError;

use crate::neighbourhood::{self, Neighbourhood};
use crate::numeric::{f64_to_f32, f64_to_u16, f64_to_usize};
use crate::sampling::fibonacci_sphere;
use std::collections::BTreeMap;

#[path = "sasa/contact.rs"]
mod contact;
pub use contact::{
    AtomContactArea, ExcludedSurfacePoint, atom_contact_areas, surface_points_excluding_pairs,
};

/// Why a solvent-accessible surface could not be computed.
#[derive(Debug, thiserror::Error)]
pub enum SasaError {
    /// There is not one radius per position.
    #[error("expected {positions} radii to match the positions, found {radii}")]
    LengthMismatch {
        /// How many positions were supplied.
        positions: usize,
        /// How many radii were supplied.
        radii: usize,
    },
    /// The probe radius is negative or not finite.
    #[error("the probe radius must be finite and non-negative")]
    InvalidProbe,
    /// A radius is negative or not finite.
    #[error("every radius must be finite and non-negative")]
    InvalidRadius,
    /// Surface sampling density is not finite and positive.
    #[error("surface density must be finite and positive")]
    InvalidDensity,
    /// No test points were requested, so no area could be sampled.
    #[error("at least one test point is required")]
    NoPoints,
    /// Grid resolution or allocation ceiling is invalid.
    #[error("grid resolution must be positive and finite and max_cells must be non-zero")]
    InvalidGridOptions,
    /// The requested grid resolution would need more cells than allowed, which
    /// guards against exhausting memory on a tiny resolution or a huge box.
    #[error("the grid would need {cells} cells, over the limit")]
    GridTooLarge {
        /// How many cells the requested grid would hold.
        cells: usize,
    },
    /// Grid dimensions overflow the host index domain before allocation.
    #[error("the requested grid dimensions exceed the host index domain")]
    GridDimensionsOverflow,
    /// The neighbour search rejected the input.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// Per-atom solvent-accessible surface area, in the squared unit of the radii.
///
/// `radii` are the bare atomic radii; the probe is added internally, so passing
/// a zero probe measures the van der Waals surface. `points` sets how finely
/// each sphere is sampled — more points cost more and converge on the exact
/// area. An atom with no neighbours returns its full sphere area exactly,
/// independent of `points`.
///
/// Runs in `O(atoms · points · local density)` time.
///
/// # Errors
///
/// Returns [`SasaError`] for mismatched lengths, a non-finite or negative probe
/// or radius, a zero point count, or a neighbour-search failure.
///
/// # Examples
///
/// ```
/// use pdbiox_surface::shrake_rupley;
/// use core::f64::consts::PI;
///
/// // A lone atom is fully exposed: its area is that of the expanded sphere.
/// let areas = shrake_rupley(&[[0.0, 0.0, 0.0]], &[2.0], 1.0, 200)?;
/// let expanded = 2.0 + 1.0;
/// assert!((areas[0] - 4.0 * PI * expanded * expanded).abs() < 1e-6);
/// # Ok::<(), pdbiox_surface::SasaError>(())
/// ```
pub fn shrake_rupley(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    points: u16,
) -> Result<Vec<f64>, SasaError> {
    let Some(geometry) = prepare(positions, radii, probe, points)? else {
        return Ok(Vec::new());
    };

    let per_point = 4.0 * core::f64::consts::PI / f64::from(points);

    let mut areas = Vec::with_capacity(positions.len());

    for atom in 0..positions.len() {
        areas.push(atom_area(atom, &geometry, per_point));
    }

    Ok(areas)
}

/// Computes sampled solvent-accessible area for one atom.
///
/// Atoms without neighbours take the exact `4πR²` fast path and therefore do
/// not scan sampling directions.
fn atom_area(atom: usize, geometry: &Geometry, per_point: f64) -> f64 {
    let radius = geometry.hood.expanded[atom];

    if radius <= 0.0 {
        return 0.0;
    }

    if geometry.hood.adjacency[atom].is_empty() {
        return 4.0 * core::f64::consts::PI * radius * radius;
    }

    let mut accessible = 0u32;

    geometry.for_each_exposed(atom, |_| {
        accessible += 1;
    });

    f64::from(accessible) * per_point * radius * radius
}

/// A sampled point on the molecular surface, with its outward normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfacePoint {
    /// Index of the atom the point sits on.
    pub atom: usize,
    /// The point's position.
    pub position: [f32; 3],
    /// The outward unit normal, pointing away from the atom centre.
    pub normal: [f32; 3],
}

/// Returns the accessible surface as a set of points, each with a normal.
///
/// A point is emitted for every test direction on every atom that is not buried
/// by a neighbour, so the density follows `samples`. The normal points radially
/// out from the atom the point belongs to.
///
/// Runs in `O(atoms · samples · local density)` time.
///
/// # Errors
///
/// Returns the same errors as [`shrake_rupley`].
pub fn surface_points(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    samples: u16,
) -> Result<Vec<SurfacePoint>, SasaError> {
    let Some(geometry) = prepare(positions, radii, probe, samples)? else {
        return Ok(Vec::new());
    };

    let mut points = Vec::new();

    for atom in 0..positions.len() {
        append_fixed_surface_points(&geometry, atom, &mut points);
    }

    Ok(points)
}

/// Appends all fixed-sampling exposed points for one atom.
///
/// Runtime is `O(samples · local density)` and no temporary direction vector is
/// created beyond the shared geometry directions.
fn append_fixed_surface_points(geometry: &Geometry, atom: usize, output: &mut Vec<SurfacePoint>) {
    let radius = geometry.hood.expanded[atom];

    if radius <= 0.0 {
        return;
    }

    let centre = geometry.hood.centres[atom];

    geometry.for_each_exposed(atom, |direction| {
        output.push(surface_point(atom, centre, radius, direction));
    });
}

/// Returns accessible surface points at an approximately uniform area density.
///
/// Each atom receives `ceil(4πR² × density)` deterministic Fibonacci samples,
/// where `R` includes the probe. Unlike a fixed samples-per-atom request, this
/// gives large and small atoms the same points per square ångström.
///
/// # Errors
///
/// Returns [`SasaError::InvalidDensity`] for a non-positive or non-finite
/// density, [`SasaError::GridTooLarge`] if one atom needs more than 65,535
/// samples, and the neighbourhood errors documented by [`shrake_rupley`].
pub fn surface_points_at_density(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    density: f32,
) -> Result<Vec<SurfacePoint>, SasaError> {
    validate_density(density)?;

    let Some(hood) = neighbourhood::build(positions, radii, probe)? else {
        return Ok(Vec::new());
    };

    let mut directions = BTreeMap::new();
    let mut points = Vec::new();

    for atom in 0..positions.len() {
        append_density_surface_points(&hood, atom, density, &mut directions, &mut points)?;
    }

    Ok(points)
}

/// Samples one atom at area-proportional Fibonacci density.
///
/// Directions are cached by sample count so atoms with the same area density
/// reuse the same deterministic lattice instead of recalculating it.
fn append_density_surface_points(
    hood: &Neighbourhood,
    atom: usize,
    density: f32,
    directions: &mut BTreeMap<u16, Vec<[f64; 3]>>,
    output: &mut Vec<SurfacePoint>,
) -> Result<(), SasaError> {
    let radius = hood.expanded[atom];

    if radius <= 0.0 {
        return Ok(());
    }

    let samples = samples_for_density(radius, density)?;
    let centre = hood.centres[atom];

    let directions = directions
        .entry(samples)
        .or_insert_with(|| fibonacci_sphere(samples));
    for &direction in directions.iter() {
        let point = point_on_sphere(centre, radius, direction);

        if hood.point_is_clear(point, atom) {
            output.push(SurfacePoint {
                atom,
                position: point.map(f64_to_f32),
                normal: direction.map(f64_to_f32),
            });
        }
    }

    Ok(())
}

/// Converts area density to a bounded Fibonacci sample count.
///
/// # Errors
///
/// Returns [`SasaError::GridTooLarge`] when the required count exceeds
/// `u16::MAX`.
fn samples_for_density(radius: f64, density: f32) -> Result<u16, SasaError> {
    let requested = (4.0 * core::f64::consts::PI * radius * radius * f64::from(density)).ceil();

    if requested > f64::from(u16::MAX) {
        return Err(SasaError::GridTooLarge {
            cells: f64_to_usize(requested),
        });
    }

    Ok(f64_to_u16(requested.max(1.0)))
}

/// Validates positive finite area sampling density.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_density(density: f32) -> Result<(), SasaError> {
    if density.is_finite() && density > 0.0 {
        Ok(())
    } else {
        Err(SasaError::InvalidDensity)
    }
}

/// Constructs one exposed surface point from a radial direction.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn surface_point(atom: usize, centre: [f64; 3], radius: f64, direction: [f64; 3]) -> SurfacePoint {
    let point = point_on_sphere(centre, radius, direction);

    SurfacePoint {
        atom,
        position: point.map(f64_to_f32),
        normal: direction.map(f64_to_f32),
    }
}

/// Returns one Cartesian point on a sphere.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn point_on_sphere(centre: [f64; 3], radius: f64, direction: [f64; 3]) -> [f64; 3] {
    [
        centre[0] + radius * direction[0],
        centre[1] + radius * direction[1],
        centre[2] + radius * direction[2],
    ]
}

/// The shared geometry a point-sampling surface needs: the probe-grown
/// neighbourhood plus the fixed set of sampling directions.
pub(crate) struct Geometry {
    /// Grown radii, centres and overlap adjacency, shared with other passes.
    pub(crate) hood: Neighbourhood,
    /// The near-uniform directions each sphere is dusted with.
    directions: Vec<[f64; 3]>,
}

impl Geometry {
    /// Calls `visit` with the direction of each of an atom's exposed points.
    ///
    /// Atoms without neighbours bypass point construction and occlusion tests.
    pub(crate) fn for_each_exposed(&self, atom: usize, mut visit: impl FnMut([f64; 3])) {
        let radius = self.hood.expanded[atom];

        if radius <= 0.0 {
            return;
        }

        if self.hood.adjacency[atom].is_empty() {
            for &direction in &self.directions {
                visit(direction);
            }

            return;
        }

        let centre = self.hood.centres[atom];

        for &direction in &self.directions {
            let point = point_on_sphere(centre, radius, direction);

            if self.hood.point_is_clear(point, atom) {
                visit(direction);
            }
        }
    }
}

/// Validates the input and builds the shared surface geometry.
///
/// Returns `Ok(None)` for an empty atom set — valid, but with nothing to sample.
///
/// Runtime consists of the shared neighbourhood build plus `O(samples)`
/// deterministic direction generation.
pub(crate) fn prepare(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    samples: u16,
) -> Result<Option<Geometry>, SasaError> {
    if samples == 0 {
        return Err(SasaError::NoPoints);
    }

    let Some(hood) = neighbourhood::build(positions, radii, probe)? else {
        return Ok(None);
    };

    Ok(Some(Geometry {
        hood,
        directions: fibonacci_sphere(samples),
    }))
}

#[cfg(test)]
#[path = "sasa_tests.rs"]
mod tests;
