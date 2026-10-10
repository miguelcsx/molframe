//! Lattice-basis reduction used to make the minimum-image search exact.

use crate::numeric::rounded_i64;

/// A basis whose columns are lattice vectors, in row-major `[row][column]` form.
pub(crate) type Matrix = [[f64; 3]; 3];

/// Integer counterpart of [`Matrix`].
pub(crate) type IntMatrix = [[i32; 3]; 3];

/// Upper bound on reduction sweeps. Every accepted step strictly shortens a
/// vector of a discrete lattice, so this is only a guard against rounding noise.
const MAX_SWEEPS: usize = 256;

/// A reduced basis together with the integer change of basis that produced it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub(crate) struct Reduced {
    /// Reduced basis; column `j` is `sum_k basis[:, k] * transform[k][j]`.
    pub(crate) basis: Matrix,
    /// Unimodular integer matrix relating the reduced to the original basis.
    pub(crate) transform: [[i64; 3]; 3],
}

fn column(matrix: &Matrix, index: usize) -> [f64; 3] {
    [matrix[0][index], matrix[1][index], matrix[2][index]]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

/// Size-reduces the columns against each other until no pair can be shortened,
/// then orders them by length. The result is a short, nearly orthogonal basis
/// of the same lattice, so a bounded search around the rounded fractional
/// position is cheap even for strongly skewed cells.
pub(crate) fn reduce(basis: Matrix) -> Reduced {
    let mut basis = basis;
    let mut transform = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

    for _ in 0..MAX_SWEEPS {
        let mut changed = false;
        for target in 0..3 {
            for source in 0..3 {
                if target == source {
                    continue;
                }
                let source_vector = column(&basis, source);
                let norm = dot(source_vector, source_vector);
                if norm <= 0.0 || !norm.is_finite() {
                    continue;
                }
                let ratio = dot(column(&basis, target), source_vector) / norm;
                let factor = ratio.round();
                if factor == 0.0 || !factor.is_finite() {
                    continue;
                }
                let integer = rounded_i64(factor);
                for row in 0..3 {
                    basis[row][target] -= factor * basis[row][source];
                    transform[row][target] -= integer * transform[row][source];
                }
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    Reduced { basis, transform }
}

/// Narrows a unimodular transform and its exact inverse to `i32`.
///
/// Returns `None` if an entry does not fit or the matrix is not unimodular,
/// in which case callers search in the unreduced basis.
pub(crate) fn narrow(transform: [[i64; 3]; 3]) -> Option<(IntMatrix, IntMatrix)> {
    let t = transform;
    let cofactor = |r0: usize, r1: usize, c0: usize, c1: usize| {
        t[r0][c0]
            .checked_mul(t[r1][c1])?
            .checked_sub(t[r0][c1].checked_mul(t[r1][c0])?)
    };
    let determinant = t[0][0]
        .checked_mul(cofactor(1, 2, 1, 2)?)?
        .checked_sub(t[0][1].checked_mul(cofactor(1, 2, 0, 2)?)?)?
        .checked_add(t[0][2].checked_mul(cofactor(1, 2, 0, 1)?)?)?;
    if determinant != 1 && determinant != -1 {
        return None;
    }
    // Adjugate over the determinant; `determinant` is its own reciprocal.
    let adjugate = [
        [
            cofactor(1, 2, 1, 2)?,
            -cofactor(0, 2, 1, 2)?,
            cofactor(0, 1, 1, 2)?,
        ],
        [
            -cofactor(1, 2, 0, 2)?,
            cofactor(0, 2, 0, 2)?,
            -cofactor(0, 1, 0, 2)?,
        ],
        [
            cofactor(1, 2, 0, 1)?,
            -cofactor(0, 2, 0, 1)?,
            cofactor(0, 1, 0, 1)?,
        ],
    ];
    let mut forward = [[0_i32; 3]; 3];
    let mut inverse = [[0_i32; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            forward[row][column] = i32::try_from(t[row][column]).ok()?;
            inverse[row][column] =
                i32::try_from(adjugate[row][column].checked_mul(determinant)?).ok()?;
        }
    }
    Some((forward, inverse))
}
