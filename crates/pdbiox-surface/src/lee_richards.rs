//! Solvent-accessible surface area by the Lee–Richards slice method.
//!
//! Where the Shrake–Rupley construction dusts each sphere with points, Lee and
//! Richards cut each grown sphere into thin parallel slices and, on every slice,
//! measure the arc of the atom's circle that no neighbour covers. Summing the
//! exposed arc over the slices — each weighted by the sphere's lateral area for
//! that thickness — gives the accessible area. A lone atom, occluded by nobody,
//! recovers the sphere area 4πR² exactly.
//!
//! The occlusion on a slice is the union of the angular intervals the neighbour
//! disks subtend on the atom's circle, so the cost is `O(atoms · slices · local
//! density)` with an interval merge per slice.

use core::f64::consts::PI;

use crate::accessible_area::SasaError;
use crate::neighbourhood::{self, Neighbourhood};

/// A full turn, the span of a slice circle.
const TWO_PI: f64 = 2.0 * PI;

/// Geometry of one horizontal slice through a grown atom.
#[derive(Clone, Copy)]
struct SliceCircle {
    radius: f64,
    absolute_height: f64,
}

/// Per-atom solvent-accessible surface area by the Lee–Richards slice method.
///
/// `radii` are the bare atomic radii and the probe is added internally; `slices`
/// sets how finely each sphere is cut, more slices costing more and converging
/// on the exact area. An atom with no neighbours returns its full sphere area
/// independent of `slices`.
///
/// Runs in `O(atoms · slices · local density)` time.
///
/// # Errors
///
/// Returns the same errors as the shared neighbourhood build, and
/// [`SasaError::NoPoints`] when `slices` is zero.
///
/// # Examples
///
/// ```
/// use pdbiox_surface::lee_richards;
/// use core::f64::consts::PI;
///
/// let areas = lee_richards(&[[0.0, 0.0, 0.0]], &[2.0], 1.0, 200)?;
/// let expanded = 2.0 + 1.0;
/// assert!((areas[0] - 4.0 * PI * expanded * expanded).abs() < 1e-6);
/// # Ok::<(), pdbiox_surface::SasaError>(())
/// ```
pub fn lee_richards(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    slices: u16,
) -> Result<Vec<f64>, SasaError> {
    if slices == 0 {
        return Err(SasaError::NoPoints);
    }

    let Some(hood) = neighbourhood::build(positions, radii, probe)? else {
        return Ok(Vec::new());
    };

    let mut areas = Vec::with_capacity(positions.len());

    for atom in 0..positions.len() {
        areas.push(atom_area(atom, &hood, slices));
    }

    Ok(areas)
}

/// The accessible area of one atom, integrated over its slices.
///
/// One angular-segment buffer is reused for every slice, avoiding a second
/// temporary vector in coverage merging.
fn atom_area(atom: usize, hood: &Neighbourhood, slices: u16) -> f64 {
    let radius = hood.expanded[atom];

    if radius <= 0.0 {
        return 0.0;
    }

    if hood.adjacency[atom].is_empty() {
        return sphere_area(radius);
    }

    let centre = hood.centres[atom];
    let thickness = 2.0 * radius / f64::from(slices);

    // Archimedes: the sphere's lateral area over a slice of this thickness is the
    // same at every height, so the exposed fraction is all that varies.
    let ring_area = TWO_PI * radius * thickness;

    let mut area = 0.0;
    let mut segments = Vec::new();

    for slice in 0..slices {
        let Some(circle) = slice_circle(centre[2], radius, thickness, slice) else {
            continue;
        };

        segments.clear();

        if slice_fully_covered(atom, hood, centre, circle, &mut segments) {
            continue;
        }

        let exposed = (1.0 - covered_fraction(&mut segments)).max(0.0);

        area += exposed * ring_area;
    }

    area
}

/// Returns the exact full area of a sphere.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn sphere_area(radius: f64) -> f64 {
    4.0 * PI * radius * radius
}

/// Computes the horizontal circle represented by one slice.
///
/// Returns `None` for a degenerate zero-radius circle.
fn slice_circle(centre_z: f64, radius: f64, thickness: f64, slice: u16) -> Option<SliceCircle> {
    let height = -radius + (f64::from(slice) + 0.5) * thickness;
    let circle_squared = radius * radius - height * height;

    if circle_squared <= 0.0 {
        return None;
    }

    Some(SliceCircle {
        radius: circle_squared.sqrt(),
        absolute_height: centre_z + height,
    })
}

/// Adds neighbour-occluded angular segments for one slice.
///
/// Returns `true` as soon as one neighbour covers the entire slice circle.
fn slice_fully_covered(
    atom: usize,
    hood: &Neighbourhood,
    centre: [f64; 3],
    circle: SliceCircle,
    segments: &mut Vec<(f64, f64)>,
) -> bool {
    for &neighbour in &hood.adjacency[atom] {
        if add_neighbour_coverage(hood, neighbour as usize, centre, circle, segments) {
            return true;
        }
    }

    false
}

/// Adds the angular interval covered by one neighbouring sphere.
///
/// Returns `true` only when the neighbour covers the whole target slice circle.
fn add_neighbour_coverage(
    hood: &Neighbourhood,
    neighbour: usize,
    centre: [f64; 3],
    circle: SliceCircle,
    segments: &mut Vec<(f64, f64)>,
) -> bool {
    let neighbour_radius = hood.expanded[neighbour];
    let neighbour_centre = hood.centres[neighbour];

    let rise = neighbour_centre[2] - circle.absolute_height;

    if rise.abs() >= neighbour_radius {
        return false;
    }

    let disk_squared = neighbour_radius * neighbour_radius - rise * rise;

    if disk_squared <= 0.0 {
        return false;
    }

    let disk_radius = disk_squared.sqrt();
    let dx = neighbour_centre[0] - centre[0];
    let dy = neighbour_centre[1] - centre[1];
    let separation_squared = dx * dx + dy * dy;
    let reach = circle.radius + disk_radius;

    if separation_squared >= reach * reach {
        return false;
    }

    let separation = separation_squared.sqrt();

    if separation + circle.radius <= disk_radius {
        return true;
    }

    if separation <= f64::EPSILON {
        return false;
    }

    let cosine = ((separation_squared + circle.radius * circle.radius - disk_radius * disk_radius)
        / (2.0 * separation * circle.radius))
        .clamp(-1.0, 1.0);

    let half = cosine.acos();
    let bearing = dy.atan2(dx);

    push_wrapped_interval(segments, bearing - half, bearing + half);

    false
}

/// Normalises one angular interval to `[0, 2π)` and splits it at wraparound.
///
/// No allocation occurs except growth of `segments`.
fn push_wrapped_interval(segments: &mut Vec<(f64, f64)>, start: f64, end: f64) {
    let width = (end - start).min(TWO_PI);
    let normalized = start.rem_euclid(TWO_PI);
    let finish = normalized + width;

    if finish <= TWO_PI {
        segments.push((normalized, finish));
    } else {
        segments.push((normalized, TWO_PI));
        segments.push((0.0, finish - TWO_PI));
    }
}

/// The fraction of the full circle covered by angular intervals.
///
/// `segments` is sorted and merged in place, so the function requires no second
/// temporary interval vector. Runtime is `O(S log S)`.
fn covered_fraction(segments: &mut [(f64, f64)]) -> f64 {
    if segments.is_empty() {
        return 0.0;
    }

    segments.sort_unstable_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
    });

    let first = segments[0];
    let mut covered = 0.0;
    let mut current_start = first.0;
    let mut current_end = first.1;

    for &(start, end) in &segments[1..] {
        if start > current_end {
            covered += current_end - current_start;
            current_start = start;
            current_end = end;
        } else if end > current_end {
            current_end = end;
        }
    }

    covered += current_end - current_start;

    (covered / TWO_PI).min(1.0)
}

#[cfg(test)]
#[path = "lee_richards_tests.rs"]
mod tests;
