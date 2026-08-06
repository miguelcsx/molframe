//! Quantities summed over a whole set of positions.
//!
//! Centroid, radius of gyration, inertia tensor, principal axes and shape
//! descriptors. All of them are one pass over the positions, so the cost is
//! linear in the number of atoms and independent of how they are arranged.
//!
//! Every one of them takes masses optionally. Weighting by mass and weighting
//! every atom equally answer different questions — a centre of mass is not a
//! geometric centre — and which was computed is something a result has to be
//! able to state, so it is a parameter rather than a convention.

use crate::eigen;
use crate::measure::dot;

/// The unweighted centre of a set of positions.
///
/// Returns `None` for an empty set: the centre of nothing is not the origin.
#[must_use]
pub fn centroid(positions: &[[f32; 3]]) -> Option<[f64; 3]> {
    weighted_centre(positions, None)
}

/// The centre of mass, or the geometric centre when no masses are given.
///
/// Masses shorter than the positions are an inconsistent request and yield
/// `None` rather than a centre computed over part of the set.
#[must_use]
pub fn centre_of_mass(positions: &[[f32; 3]], masses: &[f64]) -> Option<[f64; 3]> {
    weighted_centre(positions, Some(masses))
}

fn weighted_centre(positions: &[[f32; 3]], masses: Option<&[f64]>) -> Option<[f64; 3]> {
    if positions.is_empty() {
        return None;
    }
    if let Some(masses) = masses
        && masses.len() < positions.len()
    {
        return None;
    }

    let mut sum = [0.0f64; 3];
    let mut total = 0.0f64;
    for (index, position) in positions.iter().enumerate() {
        let weight = match masses {
            Some(masses) => match masses.get(index) {
                Some(mass) => *mass,
                None => return None,
            },
            None => 1.0,
        };
        for axis in 0..3 {
            sum[axis] += f64::from(position[axis]) * weight;
        }
        total += weight;
    }
    if total <= 0.0 {
        return None;
    }
    Some([sum[0] / total, sum[1] / total, sum[2] / total])
}

/// The radius of gyration about the centre.
///
/// Returns `None` for an empty set.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::radius_of_gyration;
///
/// let pair = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
/// assert!(radius_of_gyration(&pair, &[]).is_some_and(|rg| (rg - 1.0).abs() < 1e-12));
/// ```
#[must_use]
pub fn radius_of_gyration(positions: &[[f32; 3]], masses: &[f64]) -> Option<f64> {
    let use_masses = (!masses.is_empty()).then_some(masses);
    let centre = weighted_centre(positions, use_masses)?;

    let mut sum = 0.0f64;
    let mut total = 0.0f64;
    for (index, position) in positions.iter().enumerate() {
        let weight = match use_masses {
            Some(masses) => match masses.get(index) {
                Some(mass) => *mass,
                None => return None,
            },
            None => 1.0,
        };
        let offset = offset_from(*position, centre);
        sum += weight * dot(offset, offset);
        total += weight;
    }
    if total <= 0.0 {
        return None;
    }
    Some((sum / total).sqrt())
}

/// The inertia tensor about the centre.
///
/// Symmetric by construction, which is what lets the eigensolver read only its
/// upper triangle.
#[must_use]
pub fn inertia_tensor(positions: &[[f32; 3]], masses: &[f64]) -> Option<[[f64; 3]; 3]> {
    let use_masses = (!masses.is_empty()).then_some(masses);
    let centre = weighted_centre(positions, use_masses)?;

    let mut tensor = [[0.0f64; 3]; 3];
    for (index, position) in positions.iter().enumerate() {
        let weight = match use_masses {
            Some(masses) => match masses.get(index) {
                Some(mass) => *mass,
                None => return None,
            },
            None => 1.0,
        };
        let r = offset_from(*position, centre);
        let squared = dot(r, r);
        for row in 0..3 {
            for column in 0..3 {
                let kronecker = if row == column { squared } else { 0.0 };
                let Some(slot) = tensor.get_mut(row).and_then(|row| row.get_mut(column)) else {
                    continue;
                };
                *slot += weight * (kronecker - r[row] * r[column]);
            }
        }
    }
    Some(tensor)
}

/// The three principal axes, longest first, with their moments.
///
/// The axes are the inertia tensor's eigenvectors; ordering them by moment is
/// what makes "the long axis" a definite thing rather than whichever the solver
/// happened to produce first.
#[must_use]
pub fn principal_axes(positions: &[[f32; 3]], masses: &[f64]) -> Option<eigen::Decomposition<3>> {
    inertia_tensor(positions, masses).map(eigen::symmetric)
}

/// How far the shape departs from a sphere, on `[0, 1]`.
///
/// Zero is spherical. One is a straight line. Computed from the gyration
/// tensor's eigenvalues, so it says nothing about handedness — only elongation.
#[must_use]
pub fn asphericity(positions: &[[f32; 3]]) -> Option<f64> {
    let decomposition = gyration_axes(positions)?;
    let [a, b, c] = decomposition.values;
    let trace = a + b + c;
    if trace <= f64::EPSILON {
        return Some(0.0);
    }
    // The standard combination: the largest eigenvalue against the mean of the
    // other two, normalised by the trace so the result is scale-free.
    Some((a - 0.5 * (b + c)) / trace)
}

/// The gyration tensor's eigen-decomposition, largest extent first.
#[must_use]
pub fn gyration_axes(positions: &[[f32; 3]]) -> Option<eigen::Decomposition<3>> {
    let centre = weighted_centre(positions, None)?;
    let count = positions.len() as f64;

    let mut tensor = [[0.0f64; 3]; 3];
    for position in positions {
        let r = offset_from(*position, centre);
        for row in 0..3 {
            for column in 0..3 {
                let Some(slot) = tensor.get_mut(row).and_then(|row| row.get_mut(column)) else {
                    continue;
                };
                *slot += r[row] * r[column] / count;
            }
        }
    }
    Some(eigen::symmetric(tensor))
}

/// The offset of a position from a centre.
fn offset_from(position: [f32; 3], centre: [f64; 3]) -> [f64; 3] {
    [
        f64::from(position[0]) - centre[0],
        f64::from(position[1]) - centre[1],
        f64::from(position[2]) - centre[2],
    ]
}

#[cfg(test)]
#[path = "moments_tests.rs"]
mod tests;
