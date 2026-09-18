//! Stable fractional and Cartesian coordinate conversion for triclinic cells.

use molframe_core::structure::UnitCell;
use molframe_core::{Code, Diagnostic};

/// Precomputed matrices converting one unit cell's coordinate systems.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellTransform {
    fractional_to_cartesian: [[f64; 3]; 3],
    cartesian_to_fractional: [[f64; 3]; 3],
}

impl CellTransform {
    /// Builds conversion matrices for a non-degenerate cell.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when lengths, angles or cell volume are not finite
    /// and strictly valid.
    pub fn new(cell: &UnitCell) -> Result<Self, Diagnostic> {
        let [a, b, c] = cell.lengths;
        let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);
        if ![a, b, c, alpha, beta, gamma]
            .iter()
            .all(|value| value.is_finite())
            || a <= 0.0
            || b <= 0.0
            || c <= 0.0
            || !(0.0..std::f64::consts::PI).contains(&alpha)
            || !(0.0..std::f64::consts::PI).contains(&beta)
            || !(0.0..std::f64::consts::PI).contains(&gamma)
        {
            return Err(Diagnostic::new(Code::E5004)
                .with_message("unit cell is degenerate and cannot define fractional coordinates"));
        }
        let (cos_alpha, cos_beta, cos_gamma) = (alpha.cos(), beta.cos(), gamma.cos());
        let sin_gamma = gamma.sin();
        let c_y = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let c_z_squared = c.mul_add(c, -(c * cos_beta).powi(2) - c_y.powi(2));
        if sin_gamma.abs() <= f64::EPSILON || c_z_squared <= 0.0 {
            return Err(
                Diagnostic::new(Code::E5004).with_message("unit cell has zero or imaginary volume")
            );
        }
        let forward = [
            [a, b * cos_gamma, c * cos_beta],
            [0.0, b * sin_gamma, c_y],
            [0.0, 0.0, c_z_squared.sqrt()],
        ];
        let inverse = inverse_upper_triangular(&forward).ok_or_else(|| {
            Diagnostic::new(Code::E5004).with_message("unit cell conversion matrix is singular")
        })?;
        Ok(Self {
            fractional_to_cartesian: forward,
            cartesian_to_fractional: inverse,
        })
    }

    /// Converts a fractional position to Cartesian ångström coordinates.
    #[must_use]
    pub fn to_cartesian(&self, fractional: [f64; 3]) -> [f64; 3] {
        multiply(&self.fractional_to_cartesian, fractional)
    }

    /// Converts Cartesian ångström coordinates to fractional coordinates.
    #[must_use]
    pub fn to_fractional(&self, cartesian: [f64; 3]) -> [f64; 3] {
        multiply(&self.cartesian_to_fractional, cartesian)
    }

    /// Fractional-to-Cartesian matrix, stored by rows.
    #[must_use]
    pub const fn forward_matrix(&self) -> &[[f64; 3]; 3] {
        &self.fractional_to_cartesian
    }

    /// Cartesian-to-fractional matrix, stored by rows.
    #[must_use]
    pub const fn inverse_matrix(&self) -> &[[f64; 3]; 3] {
        &self.cartesian_to_fractional
    }
}

fn multiply(matrix: &[[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

fn inverse_upper_triangular(matrix: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let (a, b, c) = (matrix[0][0], matrix[1][1], matrix[2][2]);
    if a == 0.0 || b == 0.0 || c == 0.0 {
        return None;
    }
    Some([
        [
            1.0 / a,
            -matrix[0][1] / (a * b),
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * b) / (a * b * c),
        ],
        [0.0, 1.0 / b, -matrix[1][2] / (b * c)],
        [0.0, 0.0, 1.0 / c],
    ])
}

#[cfg(test)]
#[path = "cell_tests.rs"]
mod tests;
