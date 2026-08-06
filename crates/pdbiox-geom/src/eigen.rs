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
//! Cost is fixed by the dimension, not by the data: at most `SWEEPS` passes over
//! the off-diagonal entries.

/// The most sweeps a decomposition performs.
///
/// Cyclic Jacobi on a matrix this small converges well inside this; the cap
/// exists so that a pathological input costs bounded time rather than looping.
const SWEEPS: usize = 24;

/// Below this, an off-diagonal entry is treated as already zero.
const TOLERANCE: f64 = 1e-14;

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
    #[must_use]
    pub fn dominant(&self) -> [f64; N] {
        let mut vector = [0.0; N];
        for (row, slot) in vector.iter_mut().enumerate() {
            *slot = at(&self.vectors, row, 0);
        }
        vector
    }

    /// The eigenvector at `position`, counting from the largest eigenvalue.
    #[must_use]
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
#[must_use]
pub fn symmetric<const N: usize>(matrix: [[f64; N]; N]) -> Decomposition<N> {
    let mut work = matrix;
    let mut vectors = identity::<N>();

    for _ in 0..SWEEPS {
        if off_diagonal_norm(&work) <= TOLERANCE {
            break;
        }
        for p in 0..N {
            for q in (p + 1)..N {
                rotate(&mut work, &mut vectors, p, q);
            }
        }
    }

    let mut values = [0.0; N];
    for (index, slot) in values.iter_mut().enumerate() {
        *slot = at(&work, index, index);
    }
    sort_descending(&mut values, &mut vectors);
    Decomposition { values, vectors }
}

/// Zeroes the `(p, q)` entry with a rotation, accumulating it into `vectors`.
fn rotate<const N: usize>(
    work: &mut [[f64; N]; N],
    vectors: &mut [[f64; N]; N],
    p: usize,
    q: usize,
) {
    let apq = at(work, p, q);
    if apq.abs() <= TOLERANCE {
        return;
    }
    let app = at(work, p, p);
    let aqq = at(work, q, q);

    // The rotation angle that annihilates the entry, computed through the
    // stable form rather than through an arctangent.
    let theta = (aqq - app) / (2.0 * apq);
    let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
    let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
    let cos = 1.0 / (t * t + 1.0).sqrt();
    let sin = t * cos;

    for k in 0..N {
        let akp = at(work, k, p);
        let akq = at(work, k, q);
        set(work, k, p, cos * akp - sin * akq);
        set(work, k, q, sin * akp + cos * akq);
    }
    for k in 0..N {
        let apk = at(work, p, k);
        let aqk = at(work, q, k);
        set(work, p, k, cos * apk - sin * aqk);
        set(work, q, k, sin * apk + cos * aqk);
    }
    for k in 0..N {
        let vkp = at(vectors, k, p);
        let vkq = at(vectors, k, q);
        set(vectors, k, p, cos * vkp - sin * vkq);
        set(vectors, k, q, sin * vkp + cos * vkq);
    }
}

/// Orders eigenvalues largest first, carrying their vectors with them.
///
/// Insertion sort: the dimension is three or four, so this is a handful of
/// comparisons and no allocation.
fn sort_descending<const N: usize>(values: &mut [f64; N], vectors: &mut [[f64; N]; N]) {
    for i in 1..N {
        let mut j = i;
        while j > 0 && at_1d(values, j - 1) < at_1d(values, j) {
            values.swap(j - 1, j);
            for row in vectors.iter_mut() {
                row.swap(j - 1, j);
            }
            j -= 1;
        }
    }
}

fn identity<const N: usize>() -> [[f64; N]; N] {
    let mut matrix = [[0.0; N]; N];
    for (index, row) in matrix.iter_mut().enumerate() {
        if let Some(slot) = row.get_mut(index) {
            *slot = 1.0;
        }
    }
    matrix
}

fn off_diagonal_norm<const N: usize>(matrix: &[[f64; N]; N]) -> f64 {
    let mut total = 0.0;
    for p in 0..N {
        for q in (p + 1)..N {
            let value = at(matrix, p, q);
            total += value * value;
        }
    }
    total.sqrt()
}

/// One entry, reading a position outside the matrix as zero.
///
/// Every index here is derived from the dimension and cannot be out of range;
/// answering zero rather than stopping keeps the arithmetic total.
fn at<const N: usize>(matrix: &[[f64; N]; N], row: usize, column: usize) -> f64 {
    match matrix.get(row).and_then(|row| row.get(column)) {
        Some(value) => *value,
        None => 0.0,
    }
}

fn at_1d<const N: usize>(values: &[f64; N], index: usize) -> f64 {
    match values.get(index) {
        Some(value) => *value,
        None => 0.0,
    }
}

fn set<const N: usize>(matrix: &mut [[f64; N]; N], row: usize, column: usize, value: f64) {
    if let Some(slot) = matrix.get_mut(row).and_then(|row| row.get_mut(column)) {
        *slot = value;
    }
}

#[cfg(test)]
#[path = "eigen_tests.rs"]
mod tests;
