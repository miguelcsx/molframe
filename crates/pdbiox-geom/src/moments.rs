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
/// Runs in `O(n)` time and `O(1)` auxiliary space.
#[must_use]
pub fn centroid(positions: &[[f32; 3]]) -> Option<[f64; 3]> {
    weighted_centre(positions, None)
}

/// The centre of mass, or the geometric centre when no masses are given.
///
/// Masses shorter than the positions are an inconsistent request and yield
/// `None` rather than a centre computed over part of the set.
///
/// Runs in `O(n)` time and `O(1)` auxiliary space.
#[must_use]
pub fn centre_of_mass(positions: &[[f32; 3]], masses: &[f64]) -> Option<[f64; 3]> {
    let use_masses = if masses.is_empty() {
        None
    } else {
        Some(masses)
    };

    weighted_centre(positions, use_masses)
}

/// Computes a geometric or mass-weighted centre.
///
/// `None` selects unit weights. `Some(masses)` requires at least one mass for
/// every position. Returns `None` for empty input or a non-positive/non-finite
/// total weight. Runs in `O(n)` time and `O(1)` auxiliary space.
fn weighted_centre(positions: &[[f32; 3]], masses: Option<&[f64]>) -> Option<[f64; 3]> {
    centre_and_total(positions, masses).map(|(centre, _)| centre)
}

/// Computes a centre together with the total weight used to form it.
///
/// Returning both values lets callers reuse the denominator in a second-moment
/// pass. Runs in `O(n)` time and `O(1)` auxiliary space.
fn centre_and_total(positions: &[[f32; 3]], masses: Option<&[f64]>) -> Option<([f64; 3], f64)> {
    if positions.is_empty() {
        return None;
    }

    match masses {
        Some(masses) => weighted_centre_and_total(positions, masses),
        None => uniform_centre_and_total(positions),
    }
}

/// Computes the unweighted centre and point count.
///
/// The caller guarantees that `positions` is non-empty. Runs in `O(n)` time
/// and `O(1)` auxiliary space.
fn uniform_centre_and_total(positions: &[[f32; 3]]) -> Option<([f64; 3], f64)> {
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;

    for position in positions {
        sum_x += f64::from(position[0]);
        sum_y += f64::from(position[1]);
        sum_z += f64::from(position[2]);
    }

    let total = positions.len() as f64;
    if !total.is_finite() || total <= 0.0 {
        return None;
    }

    Some(([sum_x / total, sum_y / total, sum_z / total], total))
}

/// Computes the mass-weighted centre and total mass.
///
/// Returns `None` when `masses` is shorter than `positions` or the total mass is
/// non-positive/non-finite. Runs in `O(n)` time and `O(1)` auxiliary space.
fn weighted_centre_and_total(positions: &[[f32; 3]], masses: &[f64]) -> Option<([f64; 3], f64)> {
    if masses.len() < positions.len() {
        return None;
    }

    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    let mut total = 0.0;

    for (position, &weight) in positions.iter().zip(masses) {
        sum_x += f64::from(position[0]) * weight;
        sum_y += f64::from(position[1]) * weight;
        sum_z += f64::from(position[2]) * weight;
        total += weight;
    }

    if !total.is_finite() || total <= 0.0 {
        return None;
    }

    Some(([sum_x / total, sum_y / total, sum_z / total], total))
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
///
/// Runs in `O(n)` time using two numerically stable passes and `O(1)` space.
#[must_use]
pub fn radius_of_gyration(positions: &[[f32; 3]], masses: &[f64]) -> Option<f64> {
    let use_masses = if masses.is_empty() {
        None
    } else {
        Some(masses)
    };
    let (centre, total) = centre_and_total(positions, use_masses)?;

    let sum = match use_masses {
        Some(masses) => weighted_squared_radius(positions, masses, centre)?,
        None => uniform_squared_radius(positions, centre),
    };

    Some((sum / total).sqrt())
}

/// Accumulates the unweighted squared distance from `centre`.
///
/// Runs in `O(n)` time and `O(1)` auxiliary space.
fn uniform_squared_radius(positions: &[[f32; 3]], centre: [f64; 3]) -> f64 {
    let mut sum = 0.0;

    for &position in positions {
        let offset = offset_from(position, centre);
        sum += dot(offset, offset);
    }

    sum
}

/// Accumulates the mass-weighted squared distance from `centre`.
///
/// Returns `None` when `masses` is shorter than `positions`. Runs in `O(n)`
/// time and `O(1)` auxiliary space.
fn weighted_squared_radius(
    positions: &[[f32; 3]],
    masses: &[f64],
    centre: [f64; 3],
) -> Option<f64> {
    if masses.len() < positions.len() {
        return None;
    }

    let mut sum = 0.0;
    for (&position, &weight) in positions.iter().zip(masses) {
        let offset = offset_from(position, centre);
        sum += weight * dot(offset, offset);
    }

    Some(sum)
}

/// The inertia tensor about the centre.
///
/// Symmetric by construction, which is what lets the eigensolver read only its
/// upper triangle.
///
/// Runs in `O(n)` time using two passes and `O(1)` auxiliary space.
#[must_use]
pub fn inertia_tensor(positions: &[[f32; 3]], masses: &[f64]) -> Option<[[f64; 3]; 3]> {
    let use_masses = if masses.is_empty() {
        None
    } else {
        Some(masses)
    };
    let (centre, _) = centre_and_total(positions, use_masses)?;

    let mut tensor = [[0.0; 3]; 3];

    match use_masses {
        Some(masses) => {
            if masses.len() < positions.len() {
                return None;
            }

            for (&position, &weight) in positions.iter().zip(masses) {
                accumulate_inertia(&mut tensor, offset_from(position, centre), weight);
            }
        }
        None => {
            for &position in positions {
                accumulate_inertia(&mut tensor, offset_from(position, centre), 1.0);
            }
        }
    }

    mirror_upper_triangle(&mut tensor);
    Some(tensor)
}

/// Adds one weighted point's six independent inertia components.
///
/// Runs in `O(1)` time and writes only the upper triangle.
#[inline]
fn accumulate_inertia(tensor: &mut [[f64; 3]; 3], r: [f64; 3], weight: f64) {
    let [x, y, z] = r;
    let xx = x * x;
    let yy = y * y;
    let zz = z * z;

    tensor[0][0] += weight * (yy + zz);
    tensor[0][1] -= weight * x * y;
    tensor[0][2] -= weight * x * z;
    tensor[1][1] += weight * (xx + zz);
    tensor[1][2] -= weight * y * z;
    tensor[2][2] += weight * (xx + yy);
}

/// Copies the strict upper triangle into the strict lower triangle.
///
/// Runs in `O(1)` time for a three-by-three matrix and allocates no memory.
#[inline]
fn mirror_upper_triangle(matrix: &mut [[f64; 3]; 3]) {
    matrix[1][0] = matrix[0][1];
    matrix[2][0] = matrix[0][2];
    matrix[2][1] = matrix[1][2];
}

/// The three principal axes, longest first, with their moments.
///
/// The axes are the inertia tensor's eigenvectors; ordering them by moment is
/// what makes "the long axis" a definite thing rather than whichever the solver
/// happened to produce first.
///
/// Runs in `O(n)` time plus a fixed-cost three-by-three decomposition.
#[must_use]
pub fn principal_axes(positions: &[[f32; 3]], masses: &[f64]) -> Option<eigen::Decomposition<3>> {
    inertia_tensor(positions, masses).map(eigen::symmetric)
}

/// How far the shape departs from a sphere, on `[0, 1]`.
///
/// Zero is spherical. One is a straight line. Computed from the gyration
/// tensor's eigenvalues, so it says nothing about handedness — only elongation.
///
/// Runs in `O(n)` time plus a fixed-cost three-by-three decomposition.
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
    Some(((a - 0.5 * (b + c)) / trace).clamp(0.0, 1.0))
}

/// The gyration tensor's eigen-decomposition, largest extent first.
///
/// Runs in `O(n)` time using two passes, followed by a fixed-cost decomposition,
/// and uses `O(1)` auxiliary space.
#[must_use]
pub fn gyration_axes(positions: &[[f32; 3]]) -> Option<eigen::Decomposition<3>> {
    let (centre, count) = centre_and_total(positions, None)?;

    let mut tensor = [[0.0; 3]; 3];

    for &position in positions {
        let [x, y, z] = offset_from(position, centre);

        tensor[0][0] += x * x;
        tensor[0][1] += x * y;
        tensor[0][2] += x * z;
        tensor[1][1] += y * y;
        tensor[1][2] += y * z;
        tensor[2][2] += z * z;
    }

    let inverse_count = count.recip();

    tensor[0][0] *= inverse_count;
    tensor[0][1] *= inverse_count;
    tensor[0][2] *= inverse_count;
    tensor[1][1] *= inverse_count;
    tensor[1][2] *= inverse_count;
    tensor[2][2] *= inverse_count;

    mirror_upper_triangle(&mut tensor);
    Some(eigen::symmetric(tensor))
}

/// The offset of a position from a centre.
///
/// Runs in `O(1)` time and allocates no heap memory.
#[inline]
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
