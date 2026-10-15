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

use super::qcp;
use crate::eigen;
pub(super) use crate::numeric::exact_count;
use crate::transform::Rigid;

mod statistics;

use statistics::{FitStatistics, fit_statistics, fit_statistics_weighted};

/// Numerical controls for quaternion rigid superposition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SuperposeOptions {
    /// Relative second-invariant threshold used to classify collinearity.
    pub collinear_relative_tolerance: f64,
    /// Convergence controls for the quaternion eigendecomposition.
    pub eigen: eigen::EigenOptions,
}

impl SuperposeOptions {
    /// Collinearity threshold in the named double-precision profile.
    pub const STANDARD_COLLINEAR_RELATIVE_TOLERANCE: f64 = 1e-12;

    /// Named double-precision quaternion-fit profile.
    #[must_use]
    pub const fn standard() -> Self {
        Self {
            collinear_relative_tolerance: Self::STANDARD_COLLINEAR_RELATIVE_TOLERANCE,
            eigen: eigen::EigenOptions::standard(),
        }
    }

    fn validate(self) -> Result<Self, SuperposeError> {
        if self.collinear_relative_tolerance.is_finite()
            && self.collinear_relative_tolerance > 0.0
            && self.collinear_relative_tolerance < 1.0
        {
            Ok(self)
        } else {
            Err(SuperposeError::InvalidOptions)
        }
    }
}

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
    /// A numerical tolerance was invalid.
    InvalidOptions,
    /// Per-pair weights are absent, negative, non-finite, or sum to nothing.
    InvalidWeights,
    /// The quaternion eigensolver refused the matrix or did not converge.
    Eigen(eigen::EigenError),
    /// The point count cannot be represented exactly by the numeric kernel.
    TooManyPoints,
}

impl std::fmt::Display for SuperposeError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LengthMismatch => {
                formatter.write_str("the two sets do not have the same number of positions")
            }
            Self::TooFewPoints => {
                formatter.write_str("fewer than three positions do not fix a rotation")
            }
            Self::Degenerate => formatter.write_str("the positions are collinear"),
            Self::InvalidOptions => formatter.write_str("a numerical tolerance was invalid"),
            Self::InvalidWeights => formatter
                .write_str("per-pair weights must be finite, non-negative, and not all zero"),
            Self::Eigen(error) => write!(formatter, "the quaternion eigensolver failed: {error}"),
            Self::TooManyPoints => {
                formatter.write_str("the point count cannot be represented exactly")
            }
        }
    }
}

impl std::error::Error for SuperposeError {}

molframe_core::diagnostic_from!(SuperposeError, |error| match error {
    SuperposeError::LengthMismatch => molframe_core::Code::E5102,
    SuperposeError::TooFewPoints => molframe_core::Code::E5001,
    SuperposeError::Degenerate | SuperposeError::Eigen(_) => molframe_core::Code::E5104,
    SuperposeError::InvalidOptions | SuperposeError::InvalidWeights => molframe_core::Code::E5101,
    SuperposeError::TooManyPoints => molframe_core::Code::E1903,
});

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
/// Runs in one `O(n)` pass plus a fixed-size decomposition and allocates no heap
/// memory.
///
/// # Errors
///
/// Returns why no single best transform exists.
///
/// # Examples
///
/// ```
/// use molframe_geom::{superpose, Rigid};
///
/// let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
/// // The same shape, moved.
/// let mobile = [[5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [5.0, 1.0, 0.0]];
///
/// let fit = superpose(&mobile, &reference)?;
/// assert!(fit.rmsd < 1e-9, "a translation should fit exactly");
/// # Ok::<(), molframe_geom::SuperposeError>(())
/// ```
pub fn superpose(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
) -> Result<Superposition, SuperposeError> {
    superpose_with_options(mobile, reference, SuperposeOptions::standard())
}

/// Finds the best rigid transform with explicit numerical controls.
///
/// # Errors
///
/// Returns why no single converged transform exists.
pub fn superpose_with_options(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    options: SuperposeOptions,
) -> Result<Superposition, SuperposeError> {
    superpose_weighted_with_options(mobile, reference, &[], options)
}

/// Validates the inputs and accumulates weighted fit statistics.
fn prepared_statistics_weighted(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    weights: &[f64],
    options: SuperposeOptions,
) -> Result<(FitStatistics, f64), SuperposeError> {
    let options = options.validate()?;
    if mobile.len() != reference.len() {
        return Err(SuperposeError::LengthMismatch);
    }
    if !weights.is_empty() && weights.len() != mobile.len() {
        return Err(SuperposeError::InvalidWeights);
    }

    if mobile.len() < 3 {
        return Err(SuperposeError::TooFewPoints);
    }

    // A weight is a multiplicity: the weighted fit is the unweighted fit of a
    // set in which each point appears `w_i` times. Total weight plays the role of
    // the count, and it has to be a real, positive mass or there is nothing to
    // minimise.
    let total_weight: f64 = if weights.is_empty() {
        exact_count(mobile.len()).ok_or(SuperposeError::TooManyPoints)?
    } else {
        if weights
            .iter()
            .any(|weight| !weight.is_finite() || *weight < 0.0)
        {
            return Err(SuperposeError::InvalidWeights);
        }
        weights.iter().sum()
    };
    let count = exact_count(mobile.len()).ok_or(SuperposeError::TooManyPoints)?;
    if !total_weight.is_finite() || total_weight <= 0.0 {
        return Err(SuperposeError::InvalidWeights);
    }
    let statistics = fit_statistics_weighted(mobile, reference, weights);
    if !statistics
        .mobile_centre
        .iter()
        .all(|value| value.is_finite())
        || !statistics
            .reference_centre
            .iter()
            .all(|value| value.is_finite())
    {
        return Err(SuperposeError::InvalidWeights);
    }

    if is_collinear(
        &statistics.mobile_scatter,
        options.collinear_relative_tolerance,
    ) || is_collinear(
        &statistics.reference_scatter,
        options.collinear_relative_tolerance,
    ) {
        return Err(SuperposeError::Degenerate);
    }

    Ok((statistics, count))
}

/// Finds the rigid transform minimising the weighted sum of squared deviations.
///
/// An empty `weights` means every pair weighs one, which is exactly
/// [`superpose`]. A pair with weight zero does not enter the fit at all, so a
/// caller whose weights single out one part of the structure gets a transform
/// fitted to that part instead of a compromise over all of it.
///
/// Weights are read as non-negative and finite, and at least one must be
/// positive; anything else has no minimiser to return.
///
/// # Errors
///
/// Returns [`SuperposeError::InvalidWeights`] for weights outside that contract,
/// and otherwise the failures of [`superpose`].
pub fn superpose_weighted(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    weights: &[f64],
) -> Result<Superposition, SuperposeError> {
    superpose_weighted_with_options(mobile, reference, weights, SuperposeOptions::standard())
}

/// [`superpose_weighted`] with explicit numerical controls.
///
/// # Errors
///
/// Returns why no single converged transform exists, or why the weights cannot
/// define a minimiser.
pub fn superpose_weighted_with_options(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    weights: &[f64],
    options: SuperposeOptions,
) -> Result<Superposition, SuperposeError> {
    let (statistics, count) = prepared_statistics_weighted(mobile, reference, weights, options)?;
    let key = key_matrix(&statistics.covariance);
    let initial = initial_correlation(&statistics);
    check_eigen_options(options.eigen)?;
    // QCP first; a repeated top eigenvalue (or no convergence) falls back to
    // the general eigensolver, which `options.eigen` controls.
    let (maximum_correlation, quaternion) =
        if let Some(pair) = qcp::dominant_eigenpair(&key, initial) {
            pair
        } else {
            let decomposition =
                eigen::symmetric_with_options(key, options.eigen).map_err(SuperposeError::Eigen)?;
            (decomposition.values[0], decomposition.dominant())
        };

    let Some(rotation) = rotation_from(quaternion) else {
        return Err(SuperposeError::Degenerate);
    };

    // Rotate about the mobile centre, then move that centre onto the reference's.
    let transform = Rigid::new(
        rotation,
        translation(
            rotation,
            statistics.mobile_centre,
            statistics.reference_centre,
        ),
    );
    let rmsd = optimal_rmsd(&statistics, maximum_correlation, count);

    Ok(Superposition { transform, rmsd })
}

/// Rejects invalid eigensolver controls even when QCP never needs the solver.
fn check_eigen_options(options: eigen::EigenOptions) -> Result<(), SuperposeError> {
    options
        .validate()
        .map(|_| ())
        .map_err(SuperposeError::Eigen)
}

/// The minimum RMSD after optimal rigid superposition, without the rotation.
///
/// Only the dominant eigenvalue of the key matrix is needed, so no quaternion or
/// transform is built. Validation matches [`superpose`], and the value equals
/// `superpose(..).rmsd`. Runs in one `O(n)` pass plus a fixed-size solve and
/// allocates no heap memory.
///
/// # Errors
///
/// Returns why no single best transform exists.
pub fn rmsd_after_fit(mobile: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, SuperposeError> {
    let options = SuperposeOptions::standard();
    let (statistics, count) = prepared_statistics(mobile, reference, options)?;
    let key = key_matrix(&statistics.covariance);

    let maximum_correlation =
        if let Some(value) = qcp::max_eigenvalue(&key, initial_correlation(&statistics)) {
            value
        } else {
            eigen::symmetric_with_options(key, options.eigen)
                .map_err(SuperposeError::Eigen)?
                .values[0]
        };

    Ok(optimal_rmsd(&statistics, maximum_correlation, count))
}

/// Validates the inputs and accumulates the fit statistics.
fn prepared_statistics(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    options: SuperposeOptions,
) -> Result<(FitStatistics, f64), SuperposeError> {
    let options = options.validate()?;
    if mobile.len() != reference.len() {
        return Err(SuperposeError::LengthMismatch);
    }

    if mobile.len() < 3 {
        return Err(SuperposeError::TooFewPoints);
    }

    let count = exact_count(mobile.len()).ok_or(SuperposeError::TooManyPoints)?;
    let statistics = fit_statistics(mobile, reference);

    if is_collinear(
        &statistics.mobile_scatter,
        options.collinear_relative_tolerance,
    ) || is_collinear(
        &statistics.reference_scatter,
        options.collinear_relative_tolerance,
    ) {
        return Err(SuperposeError::Degenerate);
    }

    Ok((statistics, count))
}

/// An upper bound on the dominant eigenvalue: half the summed squared norms.
fn initial_correlation(statistics: &FitStatistics) -> f64 {
    f64::midpoint(
        trace(&statistics.mobile_scatter),
        trace(&statistics.reference_scatter),
    )
}

/// Returns whether a centred scatter tensor has rank below two.
///
/// A rank-zero or rank-one scatter represents coincident or collinear points,
/// for which rotation about one axis is unconstrained. Runs in `O(1)` time.
fn is_collinear(scatter: &[[f64; 3]; 3], relative_tolerance: f64) -> bool {
    let xx = scatter[0][0];
    let xy = scatter[0][1];
    let xz = scatter[0][2];
    let yy = scatter[1][1];
    let yz = scatter[1][2];
    let zz = scatter[2][2];
    let trace = xx + yy + zz;

    if !trace.is_finite() || trace <= 0.0 {
        return true;
    }

    let second_invariant = xx * yy + xx * zz + yy * zz - xy * xy - xz * xz - yz * yz;

    !second_invariant.is_finite() || second_invariant <= relative_tolerance * trace * trace
}

/// The symmetric four-by-four whose dominant eigenvector is the best rotation.
///
/// Runs in `O(1)` time and returns a stack-allocated matrix.
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
///
/// Returns `None` for a zero-length or non-finite quaternion. Runs in `O(1)`
/// time and allocates no memory.
fn rotation_from(quaternion: [f64; 4]) -> Option<[[f64; 3]; 3]> {
    let squared = quaternion.iter().map(|value| value * value).sum::<f64>();

    if !squared.is_finite() || squared <= 0.0 {
        return None;
    }

    let inverse_length = squared.sqrt().recip();
    let (scalar, first, second, third) = (
        quaternion[0] * inverse_length,
        quaternion[1] * inverse_length,
        quaternion[2] * inverse_length,
        quaternion[3] * inverse_length,
    );

    Some([
        [
            scalar * scalar + first * first - second * second - third * third,
            2.0 * (first * second - scalar * third),
            2.0 * (first * third + scalar * second),
        ],
        [
            2.0 * (first * second + scalar * third),
            scalar * scalar - first * first + second * second - third * third,
            2.0 * (second * third - scalar * first),
        ],
        [
            2.0 * (first * third - scalar * second),
            2.0 * (second * third + scalar * first),
            scalar * scalar - first * first - second * second + third * third,
        ],
    ])
}

/// The offset that puts the rotated mobile centre onto the reference centre.
///
/// Runs in `O(1)` time and allocates no memory.
fn translation(
    rotation: [[f64; 3]; 3],
    mobile_centre: [f64; 3],
    reference_centre: [f64; 3],
) -> [f64; 3] {
    [
        reference_centre[0]
            - (rotation[0][0] * mobile_centre[0]
                + rotation[0][1] * mobile_centre[1]
                + rotation[0][2] * mobile_centre[2]),
        reference_centre[1]
            - (rotation[1][0] * mobile_centre[0]
                + rotation[1][1] * mobile_centre[1]
                + rotation[1][2] * mobile_centre[2]),
        reference_centre[2]
            - (rotation[2][0] * mobile_centre[0]
                + rotation[2][1] * mobile_centre[1]
                + rotation[2][2] * mobile_centre[2]),
    ]
}

/// Computes the minimum RMSD from the accumulated quadratic objective.
///
/// The dominant quaternion eigenvalue is the maximum rotational correlation,
/// so no transformed point buffer or additional traversal is required. Runs in
/// `O(1)` time and allocates no memory.
fn optimal_rmsd(statistics: &FitStatistics, maximum_correlation: f64, count: f64) -> f64 {
    let mobile_squared = trace(&statistics.mobile_scatter);
    let reference_squared = trace(&statistics.reference_scatter);
    let mean_squared = (mobile_squared + reference_squared - 2.0 * maximum_correlation) / count;

    if mean_squared <= 0.0 {
        0.0
    } else {
        mean_squared.sqrt()
    }
}

/// Returns the trace of a three-by-three matrix.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn trace(matrix: &[[f64; 3]; 3]) -> f64 {
    matrix[0][0] + matrix[1][1] + matrix[2][2]
}

#[cfg(test)]
#[path = "fit_tests.rs"]
mod tests;
