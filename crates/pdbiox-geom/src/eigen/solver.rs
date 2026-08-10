//! Eigenvalues and eigenvectors of a small symmetric matrix.
//!
//! One solver serves everything here that needs one: the inertia tensor's
//! principal axes come from a three-by-three, and the rotation that best
//! superposes two point sets comes from a four-by-four. Writing it once means
//! the numerical behaviour is the same in both places.
//!
//! The method is cyclic Jacobi. For matrices this small it converges in a
//! handful of sweeps, and — the property that matters more than speed here —
//! it is deterministic: the same matrix always produces the same vectors in the
//! same order, with no dependence on iteration order or on how the values were
//! reached. A superposition whose sign flipped between runs would make every
//! coordinate downstream of it irreproducible.
//!
//! Cost is bounded by the caller-selected sweep limit.

use std::num::FpCategory;

/// Convergence controls for the cyclic-Jacobi eigensolver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EigenOptions {
    /// Relative Frobenius-norm tolerance for the off-diagonal entries.
    pub relative_tolerance: f64,
    /// Maximum complete cyclic sweeps.
    pub maximum_sweeps: usize,
}

impl EigenOptions {
    /// Relative tolerance in the named standard small-matrix profile.
    pub const STANDARD_RELATIVE_TOLERANCE: f64 = 1e-14;
    /// Sweep ceiling in the named standard small-matrix profile.
    pub const STANDARD_MAXIMUM_SWEEPS: usize = 24;

    /// Deterministic bounded profile for three- and four-dimensional matrices.
    #[must_use]
    pub const fn standard() -> Self {
        Self {
            relative_tolerance: Self::STANDARD_RELATIVE_TOLERANCE,
            maximum_sweeps: Self::STANDARD_MAXIMUM_SWEEPS,
        }
    }

    fn validate(self) -> Result<Self, EigenError> {
        if self.relative_tolerance.is_finite()
            && self.relative_tolerance > 0.0
            && self.relative_tolerance < 1.0
            && self.maximum_sweeps > 0
        {
            Ok(self)
        } else {
            Err(EigenError::InvalidOptions)
        }
    }
}

/// Why a symmetric eigendecomposition was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EigenError {
    /// The tolerance or sweep ceiling was invalid.
    InvalidOptions,
    /// At least one matrix entry was NaN or infinite.
    NonFiniteMatrix,
    /// The requested convergence was not reached before the sweep ceiling.
    DidNotConverge,
    /// An input collection was too large to count exactly in the calculation.
    InputTooLarge,
}

/// The eigenvalues and eigenvectors of a symmetric matrix.
///
/// Eigenvalues are sorted descending, and column `i` of [`Decomposition::vectors`]
/// is the unit eigenvector for eigenvalue `i`. Sorting is what makes "the largest
/// eigenvalue" and "the principal axis" well defined rather than incidental.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Decomposition<const N: usize> {
    /// Eigenvalues, descending.
    pub values: [f64; N],
    /// Eigenvectors as columns, matching `values` by position.
    pub vectors: [[f64; N]; N],
}

impl<const N: usize> Decomposition<N> {
    /// The eigenvector belonging to the largest eigenvalue.
    ///
    /// Runs in `O(N)` time and uses only the returned stack-allocated array.
    #[must_use]
    #[inline]
    pub fn dominant(&self) -> [f64; N] {
        let mut vector = [0.0; N];
        for (row, slot) in vector.iter_mut().enumerate() {
            *slot = at(&self.vectors, row, 0);
        }
        vector
    }

    /// The eigenvector at `position`, counting from the largest eigenvalue.
    ///
    /// Returns `None` when `position` is outside the decomposition. A valid
    /// lookup runs in `O(N)` time and allocates no heap memory.
    #[must_use]
    #[inline]
    pub fn vector(&self, position: usize) -> Option<[f64; N]> {
        if position >= N {
            return None;
        }

        let mut vector = [0.0; N];
        for (row, slot) in vector.iter_mut().enumerate() {
            *slot = *self.vectors.get(row)?.get(position)?;
        }
        Some(vector)
    }
}

/// Decomposes a symmetric matrix.
///
/// Only the upper triangle is read; the matrix is assumed symmetric, which every
/// caller here constructs it to be.
///
/// Uses the named [`EigenOptions::standard`] profile. Performs at most
/// `maximum_sweeps * N * (N - 1) / 2` Jacobi rotations. Space usage is
/// `O(N²)` in fixed-size stack arrays, with no heap allocation.
///
/// # Errors
///
/// Returns [`EigenError`] for non-finite input or failed convergence.
pub fn symmetric<const N: usize>(matrix: [[f64; N]; N]) -> Result<Decomposition<N>, EigenError> {
    symmetric_with_options(matrix, EigenOptions::standard())
}

/// Decomposes a symmetric matrix with explicit convergence controls.
///
/// # Errors
///
/// Returns [`EigenError`] for invalid controls, non-finite input, or failed
/// convergence.
pub fn symmetric_with_options<const N: usize>(
    matrix: [[f64; N]; N],
    options: EigenOptions,
) -> Result<Decomposition<N>, EigenError> {
    let options = options.validate()?;
    if !upper_triangle_is_finite(&matrix) {
        return Err(EigenError::NonFiniteMatrix);
    }
    let mut work = symmetric_from_upper(&matrix);
    let mut vectors = identity::<N>();
    let scale_squared = frobenius_squared_norm(&work);
    if !scale_squared.is_finite() {
        return Err(EigenError::NonFiniteMatrix);
    }
    let threshold_squared = options.relative_tolerance * options.relative_tolerance * scale_squared;
    let mut converged = off_diagonal_squared_norm(&work) <= threshold_squared;

    for _ in 0..options.maximum_sweeps {
        if converged {
            break;
        }

        for p in 0..N {
            for q in (p + 1)..N {
                rotate(&mut work, &mut vectors, p, q);
            }
        }
        converged = off_diagonal_squared_norm(&work) <= threshold_squared;
    }
    if !converged {
        return Err(EigenError::DidNotConverge);
    }

    let mut values = [0.0; N];
    for (index, slot) in values.iter_mut().enumerate() {
        *slot = at(&work, index, index);
    }

    sort_descending(&mut values, &mut vectors);
    Ok(Decomposition { values, vectors })
}

/// Zeroes the `(p, q)` entry with a rotation, accumulating it into `vectors`.
///
/// The update preserves symmetry explicitly and costs `O(N)` arithmetic with
/// no allocation. `p` and `q` are generated from the matrix dimension.
#[inline]
fn rotate<const N: usize>(
    work: &mut [[f64; N]; N],
    vectors: &mut [[f64; N]; N],
    p: usize,
    q: usize,
) {
    let apq = work[p][q];
    if matches!(apq.classify(), FpCategory::Zero) {
        return;
    }

    let app = work[p][p];
    let aqq = work[q][q];

    // The rotation angle that annihilates the entry, computed through the
    // stable form rather than through an arctangent.
    let theta = (aqq - app) / (2.0 * apq);
    let magnitude = theta.hypot(1.0);
    let t = if theta >= 0.0 {
        1.0 / (theta + magnitude)
    } else {
        -1.0 / (-theta + magnitude)
    };
    let cos = 1.0 / t.hypot(1.0);
    let sin = t * cos;

    let mut k = 0;
    while k < N {
        if k == p || k == q {
            k += 1;
            continue;
        }

        let akp = work[k][p];
        let akq = work[k][q];
        let rotated_p = cos * akp - sin * akq;
        let rotated_q = sin * akp + cos * akq;

        work[k][p] = rotated_p;
        work[p][k] = rotated_p;
        work[k][q] = rotated_q;
        work[q][k] = rotated_q;
        k += 1;
    }

    work[p][p] = app - t * apq;
    work[q][q] = aqq + t * apq;
    work[p][q] = 0.0;
    work[q][p] = 0.0;

    for row in vectors.iter_mut() {
        let vkp = row[p];
        let vkq = row[q];
        row[p] = cos * vkp - sin * vkq;
        row[q] = sin * vkp + cos * vkq;
    }
}

/// Orders eigenvalues largest first, carrying their vectors with them.
///
/// Insertion sort: the dimension is three or four, so this is a handful of
/// comparisons and no allocation.
///
/// Runs in `O(N²)` comparisons and `O(N³)` scalar swaps in the worst case;
/// for the intended dimensions this is a fixed, tiny cost.
fn sort_descending<const N: usize>(values: &mut [f64; N], vectors: &mut [[f64; N]; N]) {
    for i in 1..N {
        let mut j = i;

        while j > 0 && at_1d(values, j).total_cmp(&at_1d(values, j - 1)).is_gt() {
            values.swap(j - 1, j);
            for row in vectors.iter_mut() {
                row.swap(j - 1, j);
            }
            j -= 1;
        }
    }
}

/// Builds an `N × N` identity matrix.
///
/// Runs in `O(N²)` initialization time and uses only stack storage.
fn identity<const N: usize>() -> [[f64; N]; N] {
    let mut matrix = [[0.0; N]; N];
    for index in 0..N {
        set(&mut matrix, index, index, 1.0);
    }
    matrix
}

/// Copies the upper triangle into a complete symmetric working matrix.
///
/// The lower triangle of `matrix` is deliberately ignored. Runs in `O(N²)`
/// time and uses one fixed-size stack array.
fn symmetric_from_upper<const N: usize>(matrix: &[[f64; N]; N]) -> [[f64; N]; N] {
    let mut symmetric = [[0.0; N]; N];

    for p in 0..N {
        symmetric[p][p] = matrix[p][p];

        for q in (p + 1)..N {
            let value = matrix[p][q];
            symmetric[p][q] = value;
            symmetric[q][p] = value;
        }
    }

    symmetric
}

/// Returns the squared Frobenius norm of the strict upper triangle.
///
/// Avoiding the square root preserves the stopping condition while reducing
/// each sweep to additions and multiplications. Runs in `O(N²)` time.
fn off_diagonal_squared_norm<const N: usize>(matrix: &[[f64; N]; N]) -> f64 {
    let mut total = 0.0;

    for (p, row) in matrix.iter().enumerate() {
        for value in row.iter().skip(p + 1) {
            total += value * value;
        }
    }

    total
}

/// Squared Frobenius norm used to make convergence scale-relative.
fn frobenius_squared_norm<const N: usize>(matrix: &[[f64; N]; N]) -> f64 {
    matrix.iter().flatten().map(|value| value * value).sum()
}

fn upper_triangle_is_finite<const N: usize>(matrix: &[[f64; N]; N]) -> bool {
    matrix
        .iter()
        .enumerate()
        .all(|(row, values)| values.iter().skip(row).all(|value| value.is_finite()))
}

/// One entry, reading a position outside the matrix as zero.
///
/// Every index here is derived from the dimension and cannot be out of range;
/// answering zero rather than stopping keeps the arithmetic total.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn at<const N: usize>(matrix: &[[f64; N]; N], row: usize, column: usize) -> f64 {
    match matrix
        .get(row)
        .and_then(|matrix_row| matrix_row.get(column))
    {
        Some(value) => *value,
        None => 0.0,
    }
}

/// Reads one array entry, returning zero for an out-of-range index.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn at_1d<const N: usize>(values: &[f64; N], index: usize) -> f64 {
    match values.get(index) {
        Some(value) => *value,
        None => 0.0,
    }
}

/// Writes one matrix entry when the requested position exists.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn set<const N: usize>(matrix: &mut [[f64; N]; N], row: usize, column: usize, value: f64) {
    if let Some(slot) = matrix
        .get_mut(row)
        .and_then(|matrix_row| matrix_row.get_mut(column))
    {
        *slot = value;
    }
}

#[cfg(test)]
#[path = "solver_tests.rs"]
mod tests;
