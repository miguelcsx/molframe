//! Quaternion characteristic polynomial (QCP) solution of the key matrix.
//!
//! The key matrix `K` of a superposition is symmetric and traceless, so its
//! characteristic polynomial is the depressed quartic
//! `λ⁴ − ½·tr(K²)·λ² − ⅓·tr(K³)·λ + det K`. All four roots are real, and
//! Newton's method started at or above the largest root converges
//! monotonically to it. `(Σ|x|² + Σ|y|²) / 2` of the centred sets is such a
//! start: it bounds the largest eigenvalue from above.
//!
//! The eigenvector comes from the adjugate of `K − λI`. For a simple root that
//! matrix is rank one and proportional to `v·vᵀ`, so its largest row is the
//! eigenvector. When the top eigenvalue is (nearly) repeated the adjugate
//! vanishes and the caller must use a general eigensolver instead.

const MAX_ITERATIONS: usize = 50;
const RELATIVE_STEP: f64 = 1e-11;
const RELATIVE_GAP: f64 = 1e-10;

/// The largest eigenvalue of the key matrix, or `None` without convergence.
///
/// `initial` must be at least the largest eigenvalue. Runs in `O(1)` time and
/// allocates no memory.
pub(super) fn max_eigenvalue(k: &[[f64; 4]; 4], initial: f64) -> Option<f64> {
    let (second, third) = power_traces(k);
    let determinant = determinant(k);
    let c2 = -0.5 * second;
    let c1 = -third / 3.0;

    let mut lambda = initial;
    if !lambda.is_finite() || lambda <= 0.0 {
        return None;
    }

    for _ in 0..MAX_ITERATIONS {
        let square = lambda * lambda;
        let value = ((square + c2) * square + c1 * lambda) + determinant;
        let slope = (4.0 * square + 2.0 * c2) * lambda + c1;

        if !value.is_finite() || !slope.is_finite() || slope <= 0.0 {
            return None;
        }

        let step = value / slope;
        lambda -= step;

        if step.abs() < RELATIVE_STEP * lambda.abs() {
            return Some(lambda);
        }
    }

    None
}

/// The largest eigenvalue of the key matrix and its unit eigenvector.
///
/// Returns `None` when Newton's method does not converge or when the gap to
/// the next eigenvalue is too small for the adjugate to identify a direction.
/// Runs in `O(1)` time and allocates no memory.
pub(super) fn dominant_eigenpair(k: &[[f64; 4]; 4], initial: f64) -> Option<(f64, [f64; 4])> {
    let lambda = max_eigenvalue(k, initial)?;

    let mut shifted = *k;
    for (index, row) in shifted.iter_mut().enumerate() {
        row[index] -= lambda;
    }
    let adjugate = cofactors(&shifted);

    let scale = lambda.abs().max(frobenius(k));
    let mut best = [0.0; 4];
    let mut best_norm = 0.0;
    for row in adjugate {
        let norm = row.iter().map(|value| value * value).sum::<f64>().sqrt();
        if norm > best_norm {
            best_norm = norm;
            best = row;
        }
    }

    if !best_norm.is_finite() || best_norm < RELATIVE_GAP * scale * scale * scale {
        return None;
    }

    let inverse = best_norm.recip();
    Some((lambda, best.map(|value| value * inverse)))
}

fn frobenius(k: &[[f64; 4]; 4]) -> f64 {
    k.iter()
        .flatten()
        .map(|value| value * value)
        .sum::<f64>()
        .sqrt()
}

/// `tr(K²)` and `tr(K³)` for a symmetric matrix.
fn power_traces(k: &[[f64; 4]; 4]) -> (f64, f64) {
    let mut square = [[0.0; 4]; 4];
    for (i, row) in square.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = (0..4).map(|m| k[i][m] * k[m][j]).sum();
        }
    }
    let second = (0..4).map(|i| square[i][i]).sum();
    let third = (0..4)
        .map(|i| (0..4).map(|m| square[i][m] * k[m][i]).sum::<f64>())
        .sum();
    (second, third)
}

fn determinant(k: &[[f64; 4]; 4]) -> f64 {
    let cofactor = cofactors(k);
    (0..4).map(|j| k[0][j] * cofactor[0][j]).sum()
}

/// The cofactor matrix; its transpose is the adjugate (equal, as `K` is symmetric).
fn cofactors(a: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0; 4]; 4];
    for (i, out_row) in out.iter_mut().enumerate() {
        for (j, cell) in out_row.iter_mut().enumerate() {
            let mut minor = [[0.0; 3]; 3];
            let mut r = 0;
            for (row_index, row) in a.iter().enumerate() {
                if row_index == i {
                    continue;
                }
                let mut c = 0;
                for (column_index, value) in row.iter().enumerate() {
                    if column_index == j {
                        continue;
                    }
                    minor[r][c] = *value;
                    c += 1;
                }
                r += 1;
            }
            let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
            *cell = sign * determinant3(&minor);
        }
    }
    out
}

fn determinant3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

#[cfg(test)]
#[path = "qcp_tests.rs"]
mod tests;
