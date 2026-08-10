//! Contact patches exposed by omitting one occluding neighbour.

use super::{SasaError, point_on_sphere, samples_for_density, validate_density};
use crate::neighbourhood::{self, Neighbourhood};
use crate::numeric::f64_to_f32;
use crate::sampling::for_each_fibonacci;
use std::collections::BTreeSet;

/// A surface point exposed when one candidate neighbour is omitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExcludedSurfacePoint {
    /// Atom carrying the point.
    pub atom: usize,
    /// Candidate neighbour omitted from the accessibility test.
    pub excluded: usize,
    /// Point on the probe-grown sphere.
    pub position: [f32; 3],
    /// Surface area represented by this sample point.
    pub area: f64,
}

/// Solvent-excluded contact area for one unordered atom pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomContactArea {
    /// First atom index.
    pub first: usize,
    /// Second atom index.
    pub second: usize,
    /// Mean of the two directed buried patch areas.
    pub area: f64,
}

/// Constructs all atom-pair contact areas from probe-grown surface patches.
///
/// # Errors
///
/// Returns invalid geometry or sampling controls.
pub fn atom_contact_areas(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    density: f32,
) -> Result<Vec<AtomContactArea>, SasaError> {
    validate_density(density)?;
    let Some(hood) = neighbourhood::build(positions, radii, probe)? else {
        return Ok(Vec::new());
    };
    candidate_pairs(&hood)
        .into_iter()
        .map(|(first, second)| atom_contact_area(&hood, first, second, density))
        .collect()
}

fn candidate_pairs(hood: &Neighbourhood) -> BTreeSet<(usize, usize)> {
    let mut pairs = BTreeSet::new();
    for (first, neighbours) in hood.adjacency.iter().enumerate() {
        for &second in neighbours {
            let second = second as usize;
            if first < second {
                pairs.insert((first, second));
            }
        }
    }
    pairs
}

fn atom_contact_area(
    hood: &Neighbourhood,
    first: usize,
    second: usize,
    density: f32,
) -> Result<AtomContactArea, SasaError> {
    let mut points = Vec::new();
    sample_excluding(hood, first, second, density, &mut points)?;
    sample_excluding(hood, second, first, density, &mut points)?;
    Ok(AtomContactArea {
        first,
        second,
        area: points.iter().map(|point| point.area).sum::<f64>() * 0.5,
    })
}

/// Samples each candidate pair while omitting the partner from occlusion.
///
/// # Errors
///
/// Returns invalid geometry or sampling controls.
pub fn surface_points_excluding_pairs(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    density: f32,
    pairs: &[(usize, usize)],
) -> Result<Vec<ExcludedSurfacePoint>, SasaError> {
    validate_density(density)?;
    let Some(hood) = neighbourhood::build(positions, radii, probe)? else {
        return Ok(Vec::new());
    };
    let mut points = Vec::new();
    for &(first, second) in pairs {
        if first < positions.len() && second < positions.len() && first != second {
            sample_excluding(&hood, first, second, density, &mut points)?;
            sample_excluding(&hood, second, first, density, &mut points)?;
        }
    }
    Ok(points)
}

fn sample_excluding(
    hood: &Neighbourhood,
    atom: usize,
    excluded: usize,
    density: f32,
    output: &mut Vec<ExcludedSurfacePoint>,
) -> Result<(), SasaError> {
    let radius = hood.expanded[atom];
    if radius <= 0.0 {
        return Ok(());
    }
    let centre = hood.centres[atom];
    let samples = samples_for_density(radius, density)?;
    let area = 4.0 * core::f64::consts::PI * radius * radius / f64::from(samples);
    for_each_fibonacci(samples, |direction| {
        let point = point_on_sphere(centre, radius, direction);
        if point_inside(point, hood.centres[excluded], hood.expanded[excluded])
            && hood.point_is_clear_except(point, atom, Some(excluded))
        {
            output.push(ExcludedSurfacePoint {
                atom,
                excluded,
                position: point.map(f64_to_f32),
                area,
            });
        }
    });
    Ok(())
}

fn point_inside(point: [f64; 3], centre: [f64; 3], radius: f64) -> bool {
    point
        .iter()
        .zip(centre)
        .map(|(value, centre)| (value - centre).powi(2))
        .sum::<f64>()
        <= radius * radius
}
