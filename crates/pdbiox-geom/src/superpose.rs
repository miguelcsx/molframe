//! Fitting one set of positions onto another.
//!
//! The rotation that minimises the sum of squared deviations is found through
//! the quaternion form: build a four-by-four symmetric matrix from the
//! cross-covariance of the two centred sets, and take its dominant eigenvector.
//!
//! The quaternion route is chosen over decomposing the covariance directly for
//! one reason: it cannot produce a reflection. A decomposition can return a
//! determinant of minus one and needs a correction step, and a superposition
//! that silently mirrors a structure is exactly the kind of plausible wrong
//! answer that is hard to notice — the deviation looks fine and the chirality
//! is inverted.
//!
//! Cost is one pass over the positions plus a fixed-size decomposition, so it
//! is linear in the number of atoms.

use crate::eigen;
use crate::transform::Rigid;

/// Why a superposition could not be computed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SuperposeError {
    /// The two sets do not have the same number of positions.
    LengthMismatch,
    /// Fewer than three positions, which do not fix a rotation.
    ///
    /// Two points fix an axis and leave a free spin about it; one fixes
    /// nothing. There is no single best answer to return.
    TooFewPoints,
    /// The positions are collinear, so a rotation about their axis is free.
    Degenerate,
}

/// The deviation between two sets already in correspondence.
///
/// No fitting happens: this measures the sets where they are. Use
/// [`superpose`] first if the sets need aligning.
///
/// # Errors
///
/// Returns [`SuperposeError::LengthMismatch`] when the sets differ in size.
pub fn rmsd(mobile: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, SuperposeError> {
    if mobile.len() != reference.len() {
        return Err(SuperposeError::LengthMismatch);
    }
    if mobile.is_empty() {
        return Ok(0.0);
    }
    let mut total = 0.0f64;
    for (a, b) in mobile.iter().zip(reference) {
        total += crate::measure::distance_squared(*a, *b);
    }
    Ok((total / mobile.len() as f64).sqrt())
}

/// A fit, and what it achieved.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Superposition {
    /// The transform carrying the mobile set onto the reference.
    pub transform: Rigid,
    /// The deviation after applying it.
    pub rmsd: f64,
}

/// Finds the rigid transform that best carries `mobile` onto `reference`.
///
/// The two sets must already correspond position by position; deciding *which*
/// atom matches which is a separate question and not one geometry can answer.
///
/// # Errors
///
/// Returns why no single best transform exists.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::{superpose, Rigid};
///
/// let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
/// // The same shape, moved.
/// let mobile = [[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [5.0, 1.0, 0.0]];
///
/// let fit = superpose(&mobile, &reference)?;
/// assert!(fit.rmsd < 1e-9, "a translation should fit exactly");
/// # Ok::<(), pdbiox_geom::SuperposeError>(())
/// ```
pub fn superpose(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
) -> Result<Superposition, SuperposeError> {
    if mobile.len() != reference.len() {
        return Err(SuperposeError::LengthMismatch);
    }
    if mobile.len() < 3 {
        return Err(SuperposeError::TooFewPoints);
    }
    let Some(mobile_centre) = crate::moments::centroid(mobile) else {
        return Err(SuperposeError::TooFewPoints);
    };
    let Some(reference_centre) = crate::moments::centroid(reference) else {
        return Err(SuperposeError::TooFewPoints);
    };

    let covariance = cross_covariance(mobile, reference, mobile_centre, reference_centre);
    let decomposition = eigen::symmetric(key_matrix(&covariance));
    let quaternion = decomposition.dominant();
    let Some(rotation) = rotation_from(quaternion) else {
        return Err(SuperposeError::Degenerate);
    };

    // Rotate about the mobile centre, then move that centre onto the reference's.
    let transform = Rigid::new(
        rotation,
        translation(rotation, mobile_centre, reference_centre),
    );
    let fitted: Vec<[f32; 3]> = mobile.iter().map(|p| transform.apply(*p)).collect();
    let rmsd = rmsd(&fitted, reference)?;
    Ok(Superposition { transform, rmsd })
}

/// The three-by-three cross-covariance of the two centred sets.
fn cross_covariance(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    mobile_centre: [f64; 3],
    reference_centre: [f64; 3],
) -> [[f64; 3]; 3] {
    let mut covariance = [[0.0f64; 3]; 3];
    for (a, b) in mobile.iter().zip(reference) {
        for row in 0..3 {
            for column in 0..3 {
                let left = f64::from(a[row]) - mobile_centre[row];
                let right = f64::from(b[column]) - reference_centre[column];
                let Some(slot) = covariance.get_mut(row).and_then(|r| r.get_mut(column)) else {
                    continue;
                };
                *slot += left * right;
            }
        }
    }
    covariance
}

/// The symmetric four-by-four whose dominant eigenvector is the best rotation.
fn key_matrix(c: &[[f64; 3]; 3]) -> [[f64; 4]; 4] {
    let (xx, xy, xz) = (c[0][0], c[0][1], c[0][2]);
    let (yx, yy, yz) = (c[1][0], c[1][1], c[1][2]);
    let (zx, zy, zz) = (c[2][0], c[2][1], c[2][2]);

    [
        [xx + yy + zz, yz - zy, zx - xz, xy - yx],
        [yz - zy, xx - yy - zz, xy + yx, zx + xz],
        [zx - xz, xy + yx, yy - xx - zz, yz + zy],
        [xy - yx, zx + xz, yz + zy, zz - xx - yy],
    ]
}

/// The rotation matrix a unit quaternion describes.
fn rotation_from(q: [f64; 4]) -> Option<[[f64; 3]; 3]> {
    let length = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if length <= f64::EPSILON {
        return None;
    }
    let (w, x, y, z) = (q[0] / length, q[1] / length, q[2] / length, q[3] / length);
    Some([
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
    ])
}

/// The offset that puts the rotated mobile centre onto the reference centre.
fn translation(
    rotation: [[f64; 3]; 3],
    mobile_centre: [f64; 3],
    reference_centre: [f64; 3],
) -> [f64; 3] {
    let mut shift = [0.0f64; 3];
    for (row, slot) in shift.iter_mut().enumerate() {
        let mut rotated = 0.0;
        for column in 0..3 {
            rotated += rotation[row][column] * mobile_centre[column];
        }
        *slot = reference_centre[row] - rotated;
    }
    shift
}

#[cfg(test)]
#[path = "superpose_tests.rs"]
mod tests;
