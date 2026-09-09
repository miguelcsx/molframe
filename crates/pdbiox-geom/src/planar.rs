//! How far a set of points departs from lying in a plane.
//!
//! The best-fit plane through a set of points is the one that minimises the sum
//! of squared perpendicular distances to it, and that minimum is the smallest
//! eigenvalue of the points' scatter matrix — the same symmetric eigensolver the
//! inertia tensor uses. The root-mean-square of those distances is a single
//! number for how planar the set is: zero for a flat arrangement, growing as it
//! puckers.
//!
//! Cost is one `O(n)` pass to form the three-by-three scatter plus a fixed-size
//! decomposition.

use crate::eigen;
use crate::numeric::exact_count;

/// A plane through a set of points: a point on it and its unit normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    /// The centroid of the points, which lies on the plane.
    pub centre: [f64; 3],
    /// The unit normal to the plane.
    pub normal: [f64; 3],
}

/// Returns the RMS distance of the points from their best-fit plane.
///
/// Fewer than three points do not define a plane, so the result is `None`. A set
/// that already lies in a plane returns a value at or near zero.
///
/// Runs in `O(n)` time and `O(1)` auxiliary space.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::plane_deviation;
///
/// let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
/// let Some(deviation) = plane_deviation(&flat)? else { panic!("plane absent") };
/// assert!(deviation < 1e-6);
/// # Ok::<(), pdbiox_geom::EigenError>(())
/// ```
///
/// # Errors
///
/// Returns [`eigen::EigenError`] when the scatter matrix cannot be decomposed
/// with the named standard profile.
pub fn plane_deviation(points: &[[f32; 3]]) -> Result<Option<f64>, eigen::EigenError> {
    plane_deviation_with_options(points, eigen::EigenOptions::standard())
}

/// Returns planar deviation with explicit eigensolver convergence controls.
///
/// # Errors
///
/// Returns [`eigen::EigenError`] when the scatter matrix cannot be decomposed.
pub fn plane_deviation_with_options(
    points: &[[f32; 3]],
    options: eigen::EigenOptions,
) -> Result<Option<f64>, eigen::EigenError> {
    let Some((_, scatter, count)) = centred_scatter(points) else {
        return Ok(None);
    };
    let decomposition = eigen::symmetric_with_options(scatter, options)?;
    let smallest = decomposition.values[2];
    Ok(Some((smallest.max(0.0) / count).sqrt()))
}

/// Returns the best-fit plane through the points: their centroid and its normal.
///
/// Fewer than three points do not define a plane, so the result is `None`. The
/// normal points along the least-spread direction and is returned as a unit
/// vector; its sign is not meaningful.
///
/// Runs in `O(n)` time and `O(1)` auxiliary space.
///
/// # Errors
///
/// Returns [`eigen::EigenError`] when the scatter matrix cannot be decomposed
/// with the named standard profile.
pub fn best_fit_plane(points: &[[f32; 3]]) -> Result<Option<Plane>, eigen::EigenError> {
    best_fit_plane_with_options(points, eigen::EigenOptions::standard())
}

/// Returns a best-fit plane with explicit eigensolver convergence controls.
///
/// # Errors
///
/// Returns [`eigen::EigenError`] when the scatter matrix cannot be decomposed.
pub fn best_fit_plane_with_options(
    points: &[[f32; 3]],
    options: eigen::EigenOptions,
) -> Result<Option<Plane>, eigen::EigenError> {
    let Some((centre, scatter, _)) = centred_scatter(points) else {
        return Ok(None);
    };
    let decomposition = eigen::symmetric_with_options(scatter, options)?;
    let normal = [
        decomposition.vectors[0][2],
        decomposition.vectors[1][2],
        decomposition.vectors[2][2],
    ];
    let length = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    if length <= 0.0 {
        return Ok(None);
    }
    let inverse = length.recip();
    Ok(Some(Plane {
        centre,
        normal: [
            normal[0] * inverse,
            normal[1] * inverse,
            normal[2] * inverse,
        ],
    }))
}

/// Computes the centroid and the centred scatter matrix of a point set.
///
/// Returns `None` for fewer than three points, which do not define a plane.
fn centred_scatter(points: &[[f32; 3]]) -> Option<([f64; 3], [[f64; 3]; 3], f64)> {
    if points.len() < 3 {
        return None;
    }
    let count = exact_count(points.len())?;

    let mut mean = [0.0f64; 3];
    for point in points {
        mean[0] += f64::from(point[0]);
        mean[1] += f64::from(point[1]);
        mean[2] += f64::from(point[2]);
    }
    mean[0] /= count;
    mean[1] /= count;
    mean[2] /= count;

    let mut scatter = [[0.0f64; 3]; 3];
    for point in points {
        let x = f64::from(point[0]) - mean[0];
        let y = f64::from(point[1]) - mean[1];
        let z = f64::from(point[2]) - mean[2];
        scatter[0][0] += x * x;
        scatter[0][1] += x * y;
        scatter[0][2] += x * z;
        scatter[1][1] += y * y;
        scatter[1][2] += y * z;
        scatter[2][2] += z * z;
    }
    scatter[1][0] = scatter[0][1];
    scatter[2][0] = scatter[0][2];
    scatter[2][1] = scatter[1][2];
    Some((mean, scatter, count))
}

#[cfg(test)]
#[path = "planar_tests.rs"]
mod tests;
