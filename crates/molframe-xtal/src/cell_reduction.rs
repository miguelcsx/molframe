//! Niggli reduction of a unit cell.
//!
//! The reduced cell is the unique primitive cell of a lattice whose three edges
//! are shortest and whose angles are all acute or all obtuse. The search works on
//! the six scalar products of the edges (A, B, C, ξ, η, ζ) with the algorithms of
//! Gruber (1973) for normalization and Křivý and Gruber (1976) for reduction, with
//! the tolerance handling of Grosse-Kunstleve, Sauter and Adams (2004).

use molframe_core::structure::UnitCell;

/// Why a cell could not be reduced.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq)]
pub enum CellReductionError {
    /// A length or angle is not finite, or does not describe a lattice.
    #[error("the unit cell does not describe a lattice")]
    InvalidCell,
    /// The comparison tolerance is negative or not finite.
    #[error("the reduction tolerance must be finite and non-negative")]
    InvalidTolerance,
}

/// A cell in its reduced setting, and how to reach it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReducedCell {
    /// Edge lengths of the reduced cell in ångström.
    pub lengths: [f64; 3],
    /// Angles `α, β, γ` of the reduced cell in degrees.
    pub angles: [f64; 3],
    /// Columns are the reduced edges written in the original edges, so that
    /// `reduced = original · change_of_basis`; its determinant is one.
    pub change_of_basis: [[i32; 3]; 3],
    /// How many normalize-and-reduce rounds ran.
    pub iterations: usize,
    /// Whether the last round found nothing left to reduce.
    pub converged: bool,
}

/// Reduces a cell to its Niggli setting.
///
/// A cell that is already reduced comes back unchanged, with the identity as its
/// change of basis. `epsilon` is the absolute tolerance on the scalar products
/// (ångström squared); it decides which near-ties count as ties.
///
/// # Errors
///
/// Returns an error for a cell that does not describe a lattice, or a tolerance
/// that is negative or not finite.
pub fn niggli_reduce(
    cell: &UnitCell,
    epsilon: f64,
    iteration_limit: usize,
) -> Result<ReducedCell, CellReductionError> {
    if !epsilon.is_finite() || epsilon < 0.0 {
        return Err(CellReductionError::InvalidTolerance);
    }
    let mut vector = Gruber::of(cell)?;
    let mut iterations = 0;
    let mut converged = false;
    while iterations < iteration_limit {
        vector.normalize(epsilon);
        iterations += 1;
        // The round that reaches the limit normalizes but does not reduce.
        if iterations == iteration_limit {
            break;
        }
        if vector.reduction_step(epsilon) {
            converged = true;
            break;
        }
    }
    let (lengths, angles) = vector.cell();
    Ok(ReducedCell {
        lengths,
        angles,
        change_of_basis: vector.basis,
        iterations,
        converged,
    })
}

/// The six scalar products of a lattice basis, and the basis they describe.
#[derive(Clone, Copy, Debug)]
struct Gruber {
    a: f64,
    b: f64,
    c: f64,
    xi: f64,
    eta: f64,
    zeta: f64,
    basis: [[i32; 3]; 3],
}

impl Gruber {
    fn of(cell: &UnitCell) -> Result<Self, CellReductionError> {
        let [la, lb, lc] = cell.lengths;
        let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);
        let finite = cell
            .lengths
            .iter()
            .chain(&cell.angles)
            .all(|value| value.is_finite());
        if !finite || cell.lengths.iter().any(|length| *length <= 0.0) {
            return Err(CellReductionError::InvalidCell);
        }
        let (ca, cb, cg) = (alpha.cos(), beta.cos(), gamma.cos());
        let volume_factor = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg;
        if volume_factor <= 0.0 {
            return Err(CellReductionError::InvalidCell);
        }
        Ok(Self {
            a: la * la,
            b: lb * lb,
            c: lc * lc,
            xi: 2.0 * lb * lc * ca,
            eta: 2.0 * la * lc * cb,
            zeta: 2.0 * la * lb * cg,
            basis: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        })
    }

    fn cell(&self) -> ([f64; 3], [f64; 3]) {
        let (la, lb, lc) = (self.a.sqrt(), self.b.sqrt(), self.c.sqrt());
        let angle = |product: f64, first: f64, second: f64| {
            (product / (2.0 * first * second)).acos().to_degrees()
        };
        (
            [la, lb, lc],
            [
                angle(self.xi, lb, lc),
                angle(self.eta, la, lc),
                angle(self.zeta, la, lb),
            ],
        )
    }

    /// Orders the edges by length and makes the three products share a sign.
    fn normalize(&mut self, eps: f64) {
        self.order_first_two(eps);
        if self.b - self.c > eps
            || (self.b - self.c >= -eps && self.eta.abs() > self.zeta.abs() + eps)
        {
            std::mem::swap(&mut self.b, &mut self.c);
            std::mem::swap(&mut self.eta, &mut self.zeta);
            self.swap_columns_negated(1, 2);
            // The first two edges may now be out of order again; one more
            // pass restores them, and three swaps order any three edges.
            self.order_first_two(eps);
        }
        let positive = [self.xi, self.eta, self.zeta]
            .iter()
            .filter(|value| **value > eps)
            .count();
        let non_negative = [self.xi, self.eta, self.zeta]
            .iter()
            .filter(|value| **value >= -eps)
            .count();
        let sign = if positive == non_negative && positive % 2 == 1 {
            1.0
        } else {
            -1.0
        };
        for (value, column) in [(self.xi, 0), (self.eta, 1), (self.zeta, 2)] {
            if sign * value < -eps {
                self.negate_column(column);
            }
        }
        if positive != non_negative && positive % 2 == 1 {
            // The product that vanished is the one whose column flips the parity.
            let zeroed = match (self.zeta.abs() <= eps, self.eta.abs() <= eps) {
                (true, _) => 2,
                (false, true) => 1,
                (false, false) => 0,
            };
            self.negate_column(zeroed);
        }
        self.xi = self.xi.copysign(sign);
        self.eta = self.eta.copysign(sign);
        self.zeta = self.zeta.copysign(sign);
    }

    fn order_first_two(&mut self, eps: f64) {
        if self.a - self.b > eps
            || (self.a - self.b >= -eps && self.xi.abs() > self.eta.abs() + eps)
        {
            std::mem::swap(&mut self.a, &mut self.b);
            std::mem::swap(&mut self.xi, &mut self.eta);
            self.swap_columns_negated(0, 1);
        }
    }

    /// One reduction step; returns true when the cell is already reduced.
    fn reduction_step(&mut self, eps: f64) -> bool {
        let (a, b) = (self.a, self.b);
        let (xi, eta, zeta) = (self.xi, self.eta, self.zeta);
        let total = xi + eta + zeta + a + b;
        if xi.abs() > b + eps
            || (xi >= b - eps && 2.0 * eta < zeta - eps)
            || (xi <= -(b - eps) && zeta < -eps)
        {
            let sign = sign_of(xi);
            self.c += b - xi * sign;
            self.eta -= zeta * sign;
            self.xi -= 2.0 * b * sign;
            self.add_column(1, 2, -sign_int(xi));
        } else if eta.abs() > a + eps
            || (eta >= a - eps && 2.0 * xi < zeta - eps)
            || (eta <= -(a - eps) && zeta < -eps)
        {
            let sign = sign_of(eta);
            self.c += a - eta * sign;
            self.xi -= zeta * sign;
            self.eta -= 2.0 * a * sign;
            self.add_column(0, 2, -sign_int(eta));
        } else if zeta.abs() > a + eps
            || (zeta >= a - eps && 2.0 * xi < eta - eps)
            || (zeta <= -(a - eps) && eta < -eps)
        {
            let sign = sign_of(zeta);
            self.b += a - zeta * sign;
            self.xi -= eta * sign;
            self.zeta -= 2.0 * a * sign;
            self.add_column(0, 1, -sign_int(zeta));
        } else if total < -eps || (total <= eps && 2.0 * (a + eta) + zeta > eps) {
            self.c += a + b + xi + eta + zeta;
            self.xi += 2.0 * b + zeta;
            self.eta += 2.0 * a + zeta;
            self.add_column(0, 2, 1);
            self.add_column(1, 2, 1);
        } else {
            return true;
        }
        false
    }

    fn swap_columns_negated(&mut self, first: usize, second: usize) {
        for row in &mut self.basis {
            row.swap(first, second);
            for entry in row.iter_mut() {
                *entry = -*entry;
            }
        }
    }

    fn negate_column(&mut self, column: usize) {
        for row in &mut self.basis {
            row[column] = -row[column];
        }
    }

    fn add_column(&mut self, source: usize, target: usize, sign: i32) {
        for row in &mut self.basis {
            row[target] += sign * row[source];
        }
    }
}

fn sign_of(value: f64) -> f64 {
    if value >= 0.0 { 1.0 } else { -1.0 }
}

fn sign_int(value: f64) -> i32 {
    if value >= 0.0 { 1 } else { -1 }
}

#[cfg(test)]
#[path = "cell_reduction_tests.rs"]
mod tests;
